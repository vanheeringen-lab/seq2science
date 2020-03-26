def get_ftype(peak_caller):
    if "macs2" == peak_caller:
        if "--broad" in config["peak_caller"]["macs2"]:
            ftype = "broadPeak"
        else:
            ftype = "narrowPeak"
    if "genrich" == peak_caller:
        ftype = "narrowPeak"
    if "hmmratac" == peak_caller:
        ftype = "gappedPeak"
    return ftype


def get_genrich_replicates(wildcards):
    assembly_ish, sample_condition = "-".join(wildcards.fname.split('-')[:-1]), wildcards.fname.split('-')[-1]
    assembly = assembly_ish.split("/")[-1]

    if sample_condition in treps.index:
        return expand(f"{{dedup_dir}}/{wildcards.fname}.sambamba-queryname.bam", **config)
    else:
        return expand([f"{{dedup_dir}}/{assembly}-{replicate}.sambamba-queryname.bam"
        for replicate in treps_from_brep[(sample_condition, assembly)]], **config)


rule genrich_pileup:
    """
    Generate the pileup. We do this separately from peak-calling since these two processes have a very different
    computational footprint.
    """
    input:
        get_genrich_replicates
    output:
        bedgraphish=temp(expand("{result_dir}/genrich/{{fname}}.bdgish", **config)),
        log=temp(expand("{result_dir}/genrich/{{fname}}.log", **config))
    log:
        expand("{log_dir}/genrich_pileup/{{fname}}_pileup.log", **config)
    benchmark:
        expand("{benchmark_dir}/genrich_pileup/{{fname}}.benchmark.txt", **config)[0]
    conda:
        "../envs/genrich.yaml"
    params:
        config['peak_caller'].get('genrich', " ")  # TODO: move this to config.schema.yaml
    resources:
        mem_gb=8
    shell:
        """
        input=$(echo {input} | tr ' ' ',')
        Genrich -X -t $input -f {output.log} -k {output.bedgraphish} {params} -v > {log} 2>&1
        """


rule call_peak_genrich:
    """
    Call peaks with genrich based on the pileup.
    """
    input:
        log=expand("{result_dir}/genrich/{{fname}}.log", **config)
    output:
        narrowpeak=expand("{result_dir}/genrich/{{fname}}_peaks.narrowPeak", **config)
    log:
        expand("{log_dir}/call_peak_genrich/{{fname}}_peak.log", **config)
    benchmark:
        expand("{benchmark_dir}/call_peak_genrich/{{fname}}.benchmark.txt", **config)[0]
    conda:
        "../envs/genrich.yaml"
    params:
        config['peak_caller'].get('genrich', "")
    shell:
        "Genrich -P -f {input.log} -o {output.narrowpeak} {params} -v > {log} 2>&1"


def get_fastqc(wildcards):
    if config['layout'].get(wildcards.sample, False) == "SINGLE" or \
       config['layout'].get(wildcards.assembly, False) == "SINGLE":
        return expand("{qc_dir}/fastqc/{{sample}}_trimmed_fastqc.zip", **config)
    return sorted(expand("{qc_dir}/fastqc/{{sample}}_{fqext1}_trimmed_fastqc.zip", **config))


def get_macs2_bam(wildcards):
    if not config['macs2_keep_mates'] is True or config['layout'].get(wildcards.sample, False) == "SINGLE":
        return expand("{dedup_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config)
    return rules.keep_mates.output


rule macs2_callpeak:
    """
    Call peaks using macs2.
    Macs2 requires a genome size, which we estimate from the amount of unique kmers of the average read length.
    """
    input:
        bam=get_macs2_bam,
        fastqc=get_fastqc
    output:
        expand("{result_dir}/macs2/{{assembly}}-{{sample}}_{macs2_types}", **config)
    log:
        expand("{log_dir}/macs2_callpeak/{{assembly}}-{{sample}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/macs2_callpeak/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
    params:
        name=lambda wildcards, input: f"{wildcards.sample}" if config['layout'][wildcards.sample] == 'SINGLE' else \
                                      f"{wildcards.sample}_{config['fqext1']}",
        genome=f"{config['genome_dir']}/{{assembly}}/{{assembly}}.fa",
        macs_params=config['peak_caller'].get('macs2', ""),
        format=lambda wildcards: "BAMPE" if \
                                 (config['layout'][wildcards.sample] == "PAIRED" and "--shift" not in config['peak_caller'].get('macs2', "")) else \
                                 "BAM"
    conda:
        "../envs/macs2.yaml"
    shell:
        f"""
        # extract the kmer size, and get the effective genome size from it
        kmer_size=$(unzip -p {{input.fastqc}} {{params.name}}_trimmed_fastqc/fastqc_data.txt  | grep -P -o '(?<=Sequence length\\t).*' | grep -P -o '\d+$');
        GENSIZE=$(unique-kmers.py {{params.genome}} -k $kmer_size --quiet 2>&1 | grep -P -o '(?<=\.fa: ).*');
        echo "kmer size: $kmer_size, and effective genome size: $GENSIZE" >> {{log}}

        # call peaks
        macs2 callpeak --bdg -t {{input.bam}} --outdir {config['result_dir']}/macs2/ -n {{wildcards.assembly}}-{{wildcards.sample}} \
        {{params.macs_params}} -g $GENSIZE -f {{params.format}} >> {{log}} 2>&1
        """


rule keep_mates:
    input:
        expand("{dedup_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config)
    output:
        expand("{dedup_dir}/{{sample}}-mates-{{assembly}}.samtools-coordinate.bam", **config)
    log:
        expand("{log_dir}/keep_mates/{{assembly}}-{{sample}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/keep_mates/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
    run:
        from contextlib import redirect_stdout
        import pysam

        with open(str(log), 'w') as f:
            with redirect_stdout(f):
                paired =         1
                proper_pair =    2
                mate_unmapped =  8
                mate_reverse =   32
                first_in_pair =  64
                second_in_pair = 128

                bam_in = pysam.AlignmentFile(input[0], "rb")
                bam_out = pysam.AlignmentFile(output[0], "wb", template=bam_in)
                for line in bam_in:
                    line.qname = line.qname + ('\\1' if line.flag & first_in_pair else '\\2')
                    line.next_reference_id = 0
                    line.next_reference_start = 0
                    line.flag &= ~(
                                paired +
                                proper_pair +
                                mate_unmapped +
                                mate_reverse +
                                first_in_pair +
                                second_in_pair
                    )

                    bam_out.write(line)


rule hmmratac_genome_info:
    """
    Generate the 'genome info' that hmmratac requires for peak calling.
    https://github.com/LiuLabUB/HMMRATAC/issues/17
    
    TODO isnt this just .fa.sizes?
    """
    input:
        bam=expand("{dedup_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config)
    output:
        out=expand("{result_dir}/hmmratac/{{assembly}}-{{sample}}.genomesizes", **config),
        tmp1=temp(expand("{result_dir}/hmmratac/{{assembly}}-{{sample}}.tmp1", **config)),
        tmp2=temp(expand("{result_dir}/hmmratac/{{assembly}}-{{sample}}.tmp2", **config))
    log:
        expand("{log_dir}/hmmratac_genome_info/{{assembly}}-{{sample}}.log", **config)
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


config['hmmratac_types'] = ['.log', '.model', '_peaks.gappedPeak', '_summits.bed', '.bedgraph']

rule hmmratac:
    """
    Call 'peaks' with HMMRATAC.
    """
    input:
        genome_size=expand("{result_dir}/hmmratac/{{assembly}}-{{sample}}.genomesizes", **config),
        bam      =expand("{dedup_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config),
        bam_index=expand("{dedup_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam.bai", **config),
    output:
        expand("{result_dir}/hmmratac/{{assembly}}-{{sample}}{hmmratac_types}", **config)
    log:
        expand("{log_dir}/hmmratac/{{assembly}}-{{sample}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/hmmratac/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
    params:
        basename=lambda wildcards: expand(f"{{result_dir}}/hmmratac/{wildcards.assembly}-{wildcards.sample}", **config),
        hmmratac_params=config['peak_caller'].get('hmmratac', "")
    conda:
        "../envs/hmmratac.yaml"
    shell:
        """
        HMMRATAC --bedgraph true -o {params.basename} {params.hmmratac_params} -Xmx22G -b {input.bam} -i {input.bam_index} -g {input.genome_size} > {log} 2>&1
        """


if 'condition' in samples:
    if config['biological_replicates'] == 'idr':

        def get_idr_replicates(wildcards):
            reps = []
            for replicate in treps[(treps['assembly'] == wildcards.assembly) & (treps['condition'] == wildcards.condition)].index:
                reps.append(f"{{result_dir}}/{wildcards.peak_caller}/{wildcards.assembly}-{replicate}_peaks.{wildcards.ftype}")
            return reps

        rule idr:
            """
            Combine replicates based on the irreproducible discovery rate (IDR). Can only handle two replicates, not
            more, not less. For more than two replicates use fisher's method.
            """
            input:
                get_idr_replicates
            output:
                expand("{result_dir}/{{peak_caller}}/{{assembly}}-{{condition}}_peaks.{{ftype}}", **config),
            log:
                expand("{log_dir}/idr/{{assembly}}-{{condition}}-{{peak_caller}}-{{ftype}}.log", **config)
            benchmark:
                expand("{benchmark_dir}/idr/{{assembly}}-{{condition}}-{{peak_caller}}-{{ftype}}.benchmark.txt", **config)[0]
            params:
                lambda wildcards: "--rank 13" if wildcards.peak_caller == 'hmmratac' else ""
            conda:
                "../envs/idr.yaml"
            shell:
                """
                idr --samples {input} {params} --output-file {output} > {log} 2>&1
                """


    elif config.get('biological_replicates', "") == 'fisher':
        if 'macs2' in config['peak_caller']:

            rule macs_bdgcmp:
                """
                Prepare p-value files for rule macs_cmbreps
                """
                input:
                    treatment=expand("{result_dir}/macs2/{{assembly}}-{{sample}}_treat_pileup.bdg", **config),
                    control=  expand("{result_dir}/macs2/{{assembly}}-{{sample}}_control_lambda.bdg", **config)
                output:
                    temp(expand("{result_dir}/macs2/{{assembly}}-{{sample}}_pvalues.bdg", **config))
                log:
                    expand("{log_dir}/macs_bdgcmp/{{assembly}}-{{sample}}.log", **config)
                benchmark:
                    expand("{benchmark_dir}/macs_bdgcmp/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
                conda:
                    "../envs/macs2.yaml"
                shell:
                    """
                    macs2 bdgcmp -t {input.treatment} -c {input.control} -m ppois -o {output} > {log} 2>&1
                    """


            def get_macs_replicates(wildcards):
                return expand([f"{{result_dir}}/macs2/{wildcards.assembly}-{replicate}_pvalues.bdg"
                       for replicate in treps[(treps['assembly'] == wildcards.assembly) & (treps['condition'] == wildcards.condition)].index], **config)

            def get_macs_replicate(wildcards):
                """the original peakfile, to link if there is only 1 sample for a condition"""
                replicate = treps[(treps['assembly'] == wildcards.assembly) & (treps['condition'] == wildcards.condition)].index
                return expand(f"{{result_dir}}/macs2/{wildcards.assembly}-{replicate[0]}_peaks.{wildcards.ftype}", **config)

            rule macs_cmbreps:
                """
                Combine replicates through Fisher's method
                
                (Link original peakfile in replicate_processed if there is only 1 sample for a condition)
                """
                input:
                    bdgcmp=get_macs_replicates,
                    treatment=get_macs_replicate
                output:
                    tmpbdg=temp(expand("{result_dir}/macs2/{{assembly,.+(?<!_pvalues)}}-{{condition}}-{{ftype}}.bdg", **config)),
                    tmppeaks=temp(expand("{result_dir}/macs2/{{assembly}}-{{condition}}_peaks.temp.{{ftype}}", **config)),
                    peaks=expand("{result_dir}/macs2/{{assembly}}-{{condition}}_peaks.{{ftype}}", **config)
                log:
                    expand("{log_dir}/macs_cmbreps/{{assembly}}-{{condition}}-{{ftype}}.log", **config)
                benchmark:
                    expand("{benchmark_dir}/macs_cmbreps/{{assembly}}-{{condition}}-{{ftype}}.benchmark.txt", **config)[0]
                conda:
                    "../envs/macs2.yaml"
                params:
                    nr_reps=lambda wildcards, input: len(input.bdgcmp),
                    function="bdgpeakcall" if "--broad" not in config['peak_caller'].get('macs2', "") else "bdgbroadcall"
                shell:
                    """
                    if [ "{params.nr_reps}" == "1" ]; then
                        touch {output.tmpbdg} {output.tmppeaks}
                        mkdir -p $(dirname {output.peaks}); ln {input.treatment} {output.peaks}
                    else
                        macs2 cmbreps -i {input.bdgcmp} -o {output.tmpbdg} -m fisher > {log} 2>&1
                        macs2 {params.function} -i {output.tmpbdg} -o {output.tmppeaks} >> {log} 2>&1
                        cat {output.tmppeaks} | tail -n +2 > {output.peaks}
                    fi
                    """
