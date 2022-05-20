"""
All rules/logic related to filtering (sieving) after alignment to a genome should be here.
"""
from seq2science.util import sieve_bam


rule mark_duplicates:
    """
    Mark all duplicate reads in a bam file with picard MarkDuplicates
    """
    input:
        rules.samtools_presort.output,
    output:
        bam=expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config) if not sieve_bam(config) else temp(expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-coordinate-dupmarked.bam", **config)),
        metrics=expand("{qc_dir}/markdup/{{assembly}}-{{sample}}.samtools-coordinate.metrics.txt", **config),
    log:
        expand("{log_dir}/mark_duplicates/{{assembly}}-{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/mark_duplicates/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
    message: EXPLAIN["mark_duplicates"]
    params:
        config["markduplicates"],
    resources:
        mem_gb=8,
    conda:
        "../envs/picard.yaml"
    wildcard_constraints:
        sample=any_given("sample", "technical_replicates", "control"),
    shell:
        """
        picard MarkDuplicates TMP_DIR={resources.tmpdir} {params} INPUT={input} \
        OUTPUT={output.bam} METRICS_FILE={output.metrics} > {log} 2>&1
        """


if sieve_bam(config):

    # if doing tn5 shift, we need to re-sort afterwards!
    if config.get("tn5_shift"):
        shiftsieve = "-shifted"
    else:
        shiftsieve = ""

        ruleorder: sieve_bam > samtools_sort

    # the output of sieving depends on different preprocessing steps
    sieve_bam_output = {"final": f"{config['final_bam_dir']}/{{assembly}}-{{sample}}.samtools-coordinate{shiftsieve}.bam"}

    # if we filter on size, we make two files. One split on size, and one not.
    if config["filter_on_size"]:
        sieve_bam_output["allsizes"] = temp(
            f"{config['final_bam_dir']}/{{assembly}}-{{sample}}_allsizes.samtools-coordinate{shiftsieve}.sam"
        )

    # if we downsample to a maximum number of reads, we also need to store the results intermediately
    if config["subsample"] > -1:
        sieve_bam_output["subsample"] = temp(
            f"{config['final_bam_dir']}/{{assembly}}-{{sample}}_presubsample.samtools-coordinate{shiftsieve}.sam"
        )

    # now that we know the output sieve bam, we can mark as temp based on whether we do tn5 shift
    if config.get("tn5_shift"):
        sieve_bam_output["final"] = temp(sieve_bam_output["final"])

    # the flag to sieve on (bitwise flags)
    # https://en.wikipedia.org/wiki/SAM_(file_format)#Bitwise_flags
    sieve_flag = 0
    if config["only_primary_align"]:
        sieve_flag += 256
    if config["remove_dups"]:
        sieve_flag += 1024


    rule sieve_bam:
        """
        "Sieve" a bam.
    
        This rule does (at least) one of these:
            * filtering on minimum mapping quality
            * tn5 shift adjustment
            * remove multimappers
            * remove reads inside the blacklist
            * remove duplicates
            * filter paired-end reads on transcript length
            * subsample to have a maximum amount of reads (equal across samples)
    
        """
        input:
            bam=rules.mark_duplicates.output.bam,
            blacklist=rules.complement_blacklist.output,
            sizes=rules.get_genome_support_files.output.sizes,
        output:
            **sieve_bam_output,
        log:
            expand("{log_dir}/sieve_bam/{{assembly}}-{{sample}}.log", **config),
        benchmark:
            expand("{benchmark_dir}/sieve_bam/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
        message: EXPLAIN["sieve_bam"]
        conda:
            "../envs/samtools.yaml"
        threads: 2
        params:
            minqual=f"-q {config['min_mapping_quality']}",
            atacshift=(
                lambda wildcards, input: f" {config['rule_dir']}/../scripts/atacshift.pl /dev/stdin {input.sizes} | "
                if config["tn5_shift"]
                else ""
            ),
            blacklist=lambda wildcards, input: f"-L {input.blacklist}",
            prim_align=f"-F {sieve_flag}" if sieve_flag > 0 else "",
            sizesieve=(
                lambda wildcards, input, output: f""" tee {output.allsizes} | awk 'substr($0,1,1)=="@" || ($9>={config['min_template_length']} && $9<={config['max_template_length']}) || ($9<=-{config['min_template_length']} && $9>=-{config['max_template_length']})' | """
                if SAMPLEDICT[wildcards.sample]["layout"] == "PAIRED" and config["filter_on_size"]
                else ""
            ),
            sizesieve_touch=(
                lambda wildcards, input, output: f"touch {output.allsizes}"
                if SAMPLEDICT[wildcards.sample]["layout"] == "SINGLE" and config["filter_on_size"]
                else ""
            ),
            subsample=(
                lambda wildcards, input, output: f""" cat > {output.subsample}; nreads=$(samtools view -c {output.subsample}); if [ $nreads -gt {config['subsample']} ]; then samtools view -h -s $(echo $nreads | awk '{{print {config['subsample']}/$1}}') {output.subsample}; else samtools view -h {output.subsample}; fi | """
                if config["subsample"] > -1
                else ""
            ),
        shell:
            """
            samtools view -h {params.prim_align} {params.minqual} {params.blacklist} \
            {input.bam} | {params.atacshift} {params.sizesieve} {params.subsample}
            samtools view -b > {output.final} 2> {log}
    
            # single-end reads never output allsizes so just touch the file when filtering on size
            {params.sizesieve_touch}
            """


    rule samtools_sort_allsizes:
        """
        Sort the result of shiftsieving with the samtools sorter. This rule is identical to the
        samtools_sort rule except that the input must contain _allsizes, and the output is temporary.
        """
        input:
            rules.sieve_bam.output.final,
        output:
            temp(expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-{{sorting}}.bam",**config))
        log:
            expand("{log_dir}/samtools_sort/{{assembly}}-{{sample}}-samtools_{{sorting}}.log",**config),
        benchmark:
            expand("{benchmark_dir}/samtools_sort/{{assembly}}-{{sample}}-{{sorting}}.benchmark.txt",**config)[0]
        params:
            sort_order=lambda wildcards: "-n" if wildcards.sorting == "queryname" else "",
        wildcard_constraints:
            sample=f'({any_given("sample","technical_replicates","control")})(_allsizes)',
        threads: 2
        resources:
            mem_gb=config["bam_sort_mem"],
        conda:
            "../envs/samtools.yaml"
        shell:
            """
            # we set this trap to remove temp files when prematurely ending the rule
            trap "rm -f {resources.tmpdir}/{wildcards.assembly}-{wildcards.sample}.tmp*bam" INT;
            rm -f {resources.tmpdir}/{wildcards.assembly}-{wildcards.sample}.tmp*bam 2> {log}

            # RAM per thread in MB
            memory=$((1024*{resources.mem_gb}/{threads}))M

            samtools sort {params.sort_order} -@ {threads} -m $memory {input} -o {output} \
            -T {resources.tmpdir}/{wildcards.assembly}-{wildcards.sample}.tmp 2> {log}
            """


    rule sam2bam:
        """
        Convert a file in sam format into a bam format (binary)
        """
        input:
            "{filepath}.sam",
        output:
            temp("{filepath}.bam"),
        conda:
            "../envs/samtools.yaml"
        shell:
            """
            samtools view -b {input} > {output}
            """


rule samtools_index:
    """
    Create an index of a bam/cram file which can be used for e.g. visualization.
    """
    input:
        "{filepath}.{b}am",
    output:
        temp("{filepath}.{b}am.{b}ai"),
    params:
        config["samtools_index"],
    wildcard_constraints:
        b="b|cr"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools index {params} {input} {output}
        """


rule samtools_sort:
    """
    Sort the result of alignment or sieving with the samtools sorter.
    """
    input:
        rules.sieve_bam.output if sieve_bam(config) else rules.mark_duplicates.output.bam,
    output:
        expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-{{sorting}}.bam",**config),
    log:
        expand("{log_dir}/samtools_sort/{{assembly}}-{{sample}}-samtools_{{sorting}}.log",**config),
    benchmark:
        expand("{benchmark_dir}/samtools_sort/{{assembly}}-{{sample}}-{{sorting}}.benchmark.txt",**config)[0]
    params:
        sort_order=lambda wildcards: "-n" if wildcards.sorting == "queryname" else "",
    wildcard_constraints:
        sample=any_given("sample","technical_replicates","control")
    threads: 2
    resources:
        mem_gb=config["bam_sort_mem"],
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        # we set this trap to remove temp files when prematurely ending the rule
        trap "rm -f {resources.tmpdir}/{wildcards.assembly}-{wildcards.sample}.tmp*bam" INT;
        rm -f {resources.tmpdir}/{wildcards.assembly}-{wildcards.sample}.tmp*bam 2> {log}

        # RAM per thread in MB
        memory=$((1024*{resources.mem_gb}/{threads}))M

        samtools sort {params.sort_order} -@ {threads} -m $memory {input} -o {output} \
        -T {resources.tmpdir}/{wildcards.assembly}-{wildcards.sample}.tmp 2> {log}
        """


rule sambamba_sort:
    """
    Sort the result of alignment or sieving with the sambamba sorter.
    """
    input:
        rules.sieve_bam.output if sieve_bam(config) else rules.mark_duplicates.output.bam,
    output:
        temp(expand("{final_bam_dir}/{{assembly}}-{{sample}}.sambamba-{{sorting}}.bam", **config)),
    wildcard_constraints:
        sieve="|-sievsort",
    log:
        expand("{log_dir}/sambamba_sort/{{assembly}}-{{sample}}-sambamba_{{sorting}}.log", **config),
    message: EXPLAIN["sambamba_sort"]
    benchmark:
        expand("{benchmark_dir}/sambamba_sort/{{assembly}}-{{sample}}-{{sorting}}.benchmark.txt", **config)[0]
    params:
        sort_order=lambda wildcards: "-n" if "queryname" in wildcards.sorting else "",
        memory=f"-m {config['bam_sort_mem']}G",
    threads: 2
    conda:
        "../envs/sambamba.yaml"
    resources:
        mem_gb=config["bam_sort_mem"],
    shell:
        """
        sambamba view --nthreads {threads} -f bam  {input} -o /dev/stdout  2> {log} |
        sambamba sort --nthreads {threads} {params} /dev/stdin -o {output}  2> {log}
        """


rule bam2cram:
    """
    Convert a bam file to the more compressed cram format.
    """
    input:
        bam=rules.mark_duplicates.output.bam,
        assembly=rules.get_genome.output,
    output:
        expand("{final_bam_dir}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.cram", **config),
    message: EXPLAIN["bam2cram"]
    log:
        expand("{log_dir}/bam2cram/{{assembly}}-{{sample}}-{{sorter}}-{{sorting}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/bam2cram/{{assembly}}-{{sample}}-{{sorter}}-{{sorting}}.benchmark.txt", **config)[0]
    threads: 4
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -T {input.assembly} -C {input.bam} -@ {threads} > {output} 2> {log}
        """
