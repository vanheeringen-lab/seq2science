def get_ftype(peak_caller):
    """
    Get the filetype (narrowpeak, broadpeak, gappedpeak) for a peak caller.
    """
    if "macs2" == peak_caller:
        if "--broad" in config["peak_caller"]["macs2"]:
            ftype = "broadPeak"
        else:
            ftype = "narrowPeak"
    elif "genrich" == peak_caller:
        ftype = "narrowPeak"
    elif "hmmratac" == peak_caller:
        ftype = "gappedPeak"
    else:
        raise NotImplementedError()
    return ftype


def get_genrich_replicates(wildcards):
    assembly_ish, sample_condition = "-".join(wildcards.fname.split("-")[:-1]), wildcards.fname.split("-")[-1]
    assembly = assembly_ish.split("/")[-1]

    control = []
    if sample_condition in treps.index:
        if "control" in samples:
            control_name = treps.loc[sample_condition, "control"]
            if isinstance(control_name, str):  # ignore nan
                control = expand(f"{{final_bam_dir}}/{assembly}-{control_name}.sambamba-queryname.bam", **config)
        return {
            "control": control,
            "reps": expand(f"{{final_bam_dir}}/{wildcards.fname}.sambamba-queryname.bam", **config),
        }
    else:
        if "control" in samples:
            control_names = breps[breps["assembly"] == assembly].loc[sample_condition, "control"]
            if isinstance(control_names, str):  # single control, ignore nan
                control = expand(f"{{final_bam_dir}}/{assembly}-{control_names}.sambamba-queryname.bam", **config)
            if isinstance(control_names, pd.Series):  # multiple controls
                control = expand([f"{{final_bam_dir}}/{assembly}-{control_name}.sambamba-queryname.bam"
                                  for control_name in control_names], **config)
        return {
            "control": control,
            "reps": expand(
                [
                    f"{{final_bam_dir}}/{assembly}-{replicate}.sambamba-queryname.bam"
                    for replicate in treps_from_brep[(sample_condition, assembly)]
                ],
                **config,
            ),
        }


rule genrich_pileup:
    """
    Generate the pileup. We do this separately from peak-calling since these two
    processes have a very different computational footprint.
    """
    input:
        unpack(get_genrich_replicates),
    output:
        bedgraphish=temp(expand("{result_dir}/genrich/{{fname}}.bdgish", **config)),
        log=temp(expand("{result_dir}/genrich/{{fname}}.log", **config)),
    log:
        expand("{log_dir}/genrich_pileup/{{fname}}_pileup.log", **config),
    benchmark:
        expand("{benchmark_dir}/genrich_pileup/{{fname}}.benchmark.txt", **config)[0]
    conda:
        "../envs/genrich.yaml"
    params:
        params=config["peak_caller"].get("genrich", " "),
        control=lambda wildcards, input: f"-c {input.control}" if "control" in input else "",
    resources:
        mem_gb=8,
    shell:
        """
        input=$(echo {input.reps} | tr ' ' ',')
        Genrich -X -t $input -f {output.log} {params.control} -k {output.bedgraphish} \
        {params.params} -v > {log} 2>&1
        """


rule call_peak_genrich:
    """
    Call peaks with genrich based on the pileup.
    """
    input:
        log=expand("{result_dir}/genrich/{{fname}}.log", **config),
    output:
        narrowpeak=expand("{result_dir}/genrich/{{fname}}_peaks.narrowPeak", **config),
    log:
        expand("{log_dir}/call_peak_genrich/{{fname}}_peak.log", **config),
    benchmark:
        expand("{benchmark_dir}/call_peak_genrich/{{fname}}.benchmark.txt", **config)[0]
    message: explain_rule("call_peak_genrich")
    conda:
        "../envs/genrich.yaml"
    params:
        config["peak_caller"].get("genrich", ""),
    shell:
        """
        Genrich -P -f {input.log} -o {output.narrowpeak} {params} -v > {log} 2>&1
        """


def get_fastq_qc_file(wildcards):
    if config["trimmer"] == "trimgalore":
        if (
            sampledict.get(wildcards.sample, {}).get("layout") == "SINGLE"
            or sampledict.get(wildcards.assembly, {}).get("layout") == "SINGLE"
        ):
            return expand("{qc_dir}/fastqc/{{sample}}_trimmed_fastqc.zip", **config)
        return sorted(expand("{qc_dir}/fastqc/{{sample}}_{fqext1}_trimmed_fastqc.zip", **config))
    elif config["trimmer"] == "fastp":
        return expand("{qc_dir}/trimming/{{sample}}.fastp.json", **config)

def get_macs2_bam(wildcards):
    if not config["macs2_keep_mates"] is True or sampledict[wildcards.sample].get("layout") == "SINGLE":
        return expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config)
    return rules.keep_mates.output


def get_control_macs(wildcards):
    if not "control" in samples:
        return dict()

    control = treps.loc[wildcards.sample, "control"]
    if not isinstance(control, str) and math.isnan(control):
        return dict()

    if not config["macs2_keep_mates"] is True or sampledict[wildcards.sample].get("layout") == "SINGLE":
        return {"control": expand(f"{{final_bam_dir}}/{{{{assembly}}}}-{control}.samtools-coordinate.bam", **config)}
    return {"control": expand(f"{{final_bam_dir}}/{control}-mates-{{{{assembly}}}}.samtools-coordinate.bam", **config)}


if config["trimmer"] == "trimgalore":
    get_macs2_kmer = "kmer_size=$(unzip -p {input.fastq_qc} {params.name}_trimmed_fastqc/fastqc_data.txt  | grep -P -o '(?<=Sequence length\\t).*' | grep -P -o '\d+$')"
elif config["trimmer"] == "fastp":
    get_macs2_kmer = "kmer_size=$(jq -r .summary.after_filtering.read1_mean_length {input.fastq_qc})"


rule macs2_callpeak:
    """
    Call peaks using macs2.
    Macs2 requires a genome size, which we estimate from the amount of unique kmers of the average read length.
    """
    input:
        unpack(get_control_macs),
        bam=get_macs2_bam,
        fastq_qc=get_fastq_qc_file,
    output:
        expand("{result_dir}/macs2/{{assembly}}-{{sample}}_{macs2_types}", **config),
    log:
        expand("{log_dir}/macs2_callpeak/{{assembly}}-{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/macs2_callpeak/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
    message: explain_rule("macs2_callpeak")
    wildcard_constraints:
        sample=any_given("sample", "replicate")
    params:
        name=(
            lambda wildcards, input: wildcards.sample
            if sampledict[wildcards.sample]["layout"] == "SINGLE"
            else [f"{wildcards.sample}_{config['fqext1']}"]
        ),
        genome=f"{config['genome_dir']}/{{assembly}}/{{assembly}}.fa",
        macs_params=config["peak_caller"].get("macs2", ""),
        format=(
            lambda wildcards: "BAMPE"
            if (sampledict[wildcards.sample]["layout"] == "PAIRED" and "--shift" not in config["peak_caller"].get("macs2", ""))
            else "BAM"
        ),
        control=lambda wildcards, input: f"-c {input.control}" if "control" in input else "",
    resources:
        mem_gb=4,
    conda:
        "../envs/macs2.yaml"
    shell:
        # extract the kmer size, and get the effective genome size from it
        get_macs2_kmer +
        """
        echo "preparing to run unique-kmers.py with -k $kmer_size" >> {log}
        GENSIZE=$(unique-kmers.py {params.genome} -k $kmer_size --quiet 2>&1 | grep -P -o '(?<=\.fa: ).*')
        echo "kmer size: $kmer_size, and effective genome size: $GENSIZE" >> {log}

        # call peaks
        macs2 callpeak --bdg -t {input.bam} {params.control} --outdir {config[result_dir]}/macs2/ -n {wildcards.assembly}-{wildcards.sample} \
        {params.macs_params} -g $GENSIZE -f {params.format} >> {log} 2>&1
        """


rule keep_mates:
    """
    In-house script that, after alignment, removes the information that reads are paired.
    This can be beneficial when peak calling with macs2 when shifting + extending, since
    macs2 in this case only keeps the first in pair.
    """
    input:
        expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config),
    output:
        expand("{final_bam_dir}/{{assembly}}-{{sample}}-mates.samtools-coordinate.bam", **config),
    message: explain_rule("keep_mates")
    log:
        expand("{log_dir}/keep_mates/{{assembly}}-{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/keep_mates/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
    conda:
        "../envs/pysam.yaml"
    script:
        f"{config['rule_dir']}/../scripts/keep_mates.py"



rule hmmratac_genome_info:
    """
    Generate the 'genome info' that hmmratac requires for peak calling.
    https://github.com/LiuLabUB/HMMRATAC/issues/17

    TODO: isnt this just .fa.sizes?
    """
    input:
        bam=expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config),
    output:
        out=expand("{result_dir}/hmmratac/{{assembly}}-{{sample}}.genomesizes", **config),
        tmp1=temp(expand("{result_dir}/hmmratac/{{assembly}}-{{sample}}.tmp1", **config)),
        tmp2=temp(expand("{result_dir}/hmmratac/{{assembly}}-{{sample}}.tmp2", **config)),
    log:
        expand("{log_dir}/hmmratac_genome_info/{{assembly}}-{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/hmmratac_genome_info/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -H {input} 2>  {log} | grep SQ 2>> {log} | cut -f 2-3 2>> {log} | cut -d ':' -f 2   2>> {log} | cut -f 1        > {output.tmp1} 2> {log}
        samtools view -H {input} 2>> {log} | grep SQ 2>> {log} | cut -f 2-3 2>> {log} | cut -d ':' -f 2,3 2>> {log} | cut -d ':' -f 2 > {output.tmp2} 2> {log}
        paste {output.tmp1} {output.tmp2} > {output.out} 2> {log}
        """


config["hmmratac_types"] = [".log", ".model", "_peaks.gappedPeak", "_summits.bed", ".bedgraph"]


rule hmmratac:
    """
    Call gappedpeaks with HMMRATAC.
    """
    input:
        genome_size=expand("{result_dir}/hmmratac/{{assembly}}-{{sample}}.genomesizes", **config),
        bam=expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config),
        bam_index=expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam.bai", **config),
    output:
        expand("{result_dir}/hmmratac/{{assembly}}-{{sample}}{hmmratac_types}", **config),
    log:
        expand("{log_dir}/hmmratac/{{assembly}}-{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/hmmratac/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
    params:
        basename=lambda wildcards: expand(f"{{result_dir}}/hmmratac/{wildcards.assembly}-{wildcards.sample}", **config),
        hmmratac_params=config["peak_caller"].get("hmmratac", ""),
    conda:
        "../envs/hmmratac.yaml"
    shell:
        """
        HMMRATAC --bedgraph true -o {params.basename} {params.hmmratac_params} -Xmx22G -b {input.bam} -i {input.bam_index} -g {input.genome_size} > {log} 2>&1
        """


if "condition" in samples:
    if config["biological_replicates"] == "idr":
        ruleorder: idr > macs2_callpeak > call_peak_genrich

        def get_idr_replicates(wildcards):
            reps = []
            for replicate in treps[
                (treps["assembly"] == ori_assembly(wildcards.assembly)) & (treps["condition"] == wildcards.condition)
            ].index:
                reps.append(
                    f"{config['result_dir']}/{wildcards.peak_caller}/{wildcards.assembly}-{replicate}_peaks.{wildcards.ftype}"
                )
            return reps

        rule idr:
            """
            Combine replicates based on the irreproducible discovery rate (IDR). Can only handle two replicates.
            For more than two replicates use fisher's method.

            When combining narrowpeak files with IDR, the q-score and summit are set to -1 (means not set). However some
            downstream tools require these to be set. So we set the q-score to zero, and place the summit of the peak in
            the middle of the peak.
            """
            input:
                get_idr_replicates,
            output:
                true=expand("{result_dir}/{{peak_caller}}/{{assembly}}-{{condition}}_peaks.{{ftype}}", **config),
                temp=temp(expand("{result_dir}/{{peak_caller}}/{{assembly}}-{{condition}}_peaks.tmp.{{ftype}}", **config)),
            message: explain_rule("idr")
            log:
                expand("{log_dir}/idr/{{assembly}}-{{condition}}-{{peak_caller}}-{{ftype}}.log", **config),
            benchmark:
                expand("{benchmark_dir}/idr/{{assembly}}-{{condition}}-{{peak_caller}}-{{ftype}}.benchmark.txt", **config)[0]
            params:
                rank=lambda wildcards: "--rank 13" if wildcards.peak_caller == "hmmratac" else "",
                nr_reps=lambda wildcards, input: len(input),
            conda:
                "../envs/idr.yaml"
            shell:
                """
                if [ "{params.nr_reps}" == "1" ]; then
                    cp {input} {output.true}
                    touch {output.temp}
                else
                    idr --samples {input} {params.rank} --output-file {output.temp} > {log} 2>&1
                    if [ "{wildcards.ftype}" == "narrowPeak" ]; then
                        awk 'IFS="\t",OFS="\t" {{$9=0; $10=int(($3 - $2) / 2); print}}' {output.temp} > {output.true} 2>> {log}
                    else
                        cp {output.temp} {output.true}
                    fi
                fi
                """


    elif config.get("biological_replicates", "") == "fisher":
        if "macs2" in config["peak_caller"]:
            ruleorder: macs_cmbreps > macs2_callpeak > call_peak_genrich

            rule macs_bdgcmp:
                """
                Prepare p-value files for rule macs_cmbreps
                """
                input:
                    treatment=expand("{result_dir}/macs2/{{assembly}}-{{sample}}_treat_pileup.bdg", **config),
                    control=expand("{result_dir}/macs2/{{assembly}}-{{sample}}_control_lambda.bdg", **config),
                output:
                    temp(expand("{result_dir}/macs2/{{assembly}}-{{sample}}_qvalues.bdg", **config)),
                log:
                    expand("{log_dir}/macs_bdgcmp/{{assembly}}-{{sample}}.log", **config),
                benchmark:
                    expand("{benchmark_dir}/macs_bdgcmp/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
                resources:
                    mem_gb=4,
                conda:
                    "../envs/macs2.yaml"
                shell:
                    """
                    macs2 bdgcmp -t {input.treatment} -c {input.control} -m qpois -o {output} > {log} 2>&1
                    """


            def get_macs_replicates(wildcards):
                return expand(
                    [
                        f"{{result_dir}}/macs2/{wildcards.assembly}-{replicate}_qvalues.bdg"
                        for replicate in treps[
                            (treps["assembly"] == ori_assembly(wildcards.assembly)) & (treps["condition"] == wildcards.condition)
                        ].index
                    ],
                    **config,
                )

            def get_macs_replicate(wildcards):
                """the original peakfile, to link if there is only 1 sample for a condition"""
                replicate = treps[
                    (treps['assembly'] == ori_assembly(wildcards.assembly)) & (treps['condition'] == wildcards.condition)
                ].index
                return expand(
                    f"{{result_dir}}/macs2/{wildcards.assembly}-{replicate[0]}_peaks.{wildcards.ftype}",
                    **config,
                )

            rule macs_cmbreps:
                """
                Combine replicates through Fisher's method

                (Link original peakfile in replicate_processed if there is only 1 sample for a condition)
                """
                input:
                    bdgcmp=get_macs_replicates,
                    treatment=get_macs_replicate,
                output:
                    tmpbdg=temp(expand("{result_dir}/macs2/{{assembly,.+(?<!_qvalues)}}-{{condition}}-{{ftype}}.bdg", **config)),
                    tmppeaks=temp(expand("{result_dir}/macs2/{{assembly}}-{{condition}}_peaks.temp.{{ftype}}", **config)),
                    peaks=expand("{result_dir}/macs2/{{assembly}}-{{condition}}_peaks.{{ftype}}", **config),
                message: explain_rule("macs_cmbreps")
                log:
                    expand("{log_dir}/macs_cmbreps/{{assembly}}-{{condition}}-{{ftype}}.log", **config),
                benchmark:
                    expand("{benchmark_dir}/macs_cmbreps/{{assembly}}-{{condition}}-{{ftype}}.benchmark.txt", **config)[0]
                conda:
                    "../envs/macs2.yaml"
                params:
                    nr_reps=lambda wildcards, input: len(input.bdgcmp),
                    function="bdgpeakcall" if "--broad" not in config["peak_caller"].get("macs2", "") else "bdgbroadcall",
                    config=config["macs_cmbreps"],
                resources:
                    mem_gb=4,
                shell:
                    """
                    if [ "{params.nr_reps}" == "1" ]; then
                        touch {output.tmpbdg} {output.tmppeaks}
                        mkdir -p $(dirname {output.peaks}); cp {input.treatment} {output.peaks}
                    else
                        macs2 cmbreps -i {input.bdgcmp} -o {output.tmpbdg} -m fisher > {log} 2>&1
                        macs2 {params.function} {params.config} -i {output.tmpbdg} -o {output.tmppeaks} >> {log} 2>&1
                        cat {output.tmppeaks} | tail -n +2 > {output.peaks}
                    fi
                    """
