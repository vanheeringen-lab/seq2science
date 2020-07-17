def get_blacklist_files(wildcards):
    files = {}
    # ideally get genome is a checkpoint, however there are quite some Snakemake
    # bugs related to this. So for now we solve it like this
    # TODO: switch back to checkpoints
    if config.get("remove_blacklist") and wildcards.assembly.lower() in ["ce10", "dm3", "hg38", "hg19", "mm9", "mm10"]:
        blacklist = f"{config['genome_dir']}/{wildcards.assembly}/{wildcards.assembly}.fa"
        files["blacklist"] = blacklist

    if config.get("remove_mito"):
        sizes = f"{config['genome_dir']}/{wildcards.assembly}/{wildcards.assembly}.fa.sizes"
        files["sizes"] = sizes

    return files


rule setup_blacklist:
    """
    Combine the encode blacklist with mitochrondrial dna depending on config.
    """
    input:
        unpack(get_blacklist_files),
    output:
        temp(expand("{genome_dir}/{{assembly}}/{{assembly}}.customblacklist.bed", **config)),
    run:
        newblacklist = ""
        if config.get("remove_blacklist") and wildcards.assembly.lower() in ["ce10", "dm3", "hg38", "hg19", "mm9", "mm10"]:
            blacklist = f"{config['genome_dir']}/{wildcards.assembly}/{wildcards.assembly}.blacklist.bed"
            with open(blacklist) as file:
                newblacklist += file.read()

        if any(".fa.sizes" in inputfile for inputfile in input):
            with open(input.sizes, "r") as file:
                sizesfile = file.read().strip()
                for match in re.findall("chrM.*|chrm.*|MT.*", sizesfile):
                    chrm, size = match.split("\t")
                    newblacklist += f"{chrm}\t0\t{size}\n"

        with open(output[0], "w") as f:
            f.write(newblacklist)



rule complement_blacklist:
    """
    Take the complement of the blacklist. We need this complement to tell samtools
    to only keep reads that are in the complement of the blacklist.
    """
    input:
        blacklist=rules.setup_blacklist.output,
        sizes=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config),
    output:
        temp(expand("{genome_dir}/{{assembly}}/{{assembly}}.customblacklist_complement.bed", **config)),
    log:
        expand("{log_dir}/complement_blacklist/{{assembly}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/complement_blacklist/{{assembly}}.benchmark.txt", **config)[0]
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        sortBed -faidx {input.sizes} -i {input.blacklist} |
        complementBed -i stdin -g {input.sizes} > {output} 2> {log}
        """


rule sieve_bam:
    """
    "Sieve" a bam.

    This rule does (at least) one of these:
        * filtering on minimum mapping quality
        * tn5 shift adjustment
        * remove multimappers
        * remove reads inside the blacklist
    """
    input:
        bam=rules.samtools_presort.output,
        blacklist=rules.complement_blacklist.output,
        sizes=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config),
    output:
        temp(expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-coordinate-sieved.bam", **config)),
    log:
        expand("{log_dir}/sieve_bam/{{assembly}}-{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/sieve_bam/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
    message: explain_rule("sieve_bam")
    params:
        minqual=f"-q {config['min_mapping_quality']}",
        atacshift=(
            lambda wildcards, input: f"| {config['rule_dir']}/../scripts/atacshift.pl /dev/stdin {input.sizes}"
            if config["tn5_shift"]
            else ""
        ),
        blacklist=lambda wildcards, input: f"-L {input.blacklist}",
        prim_align=f"-F 256" if config["only_primary_align"] else "",
    conda:
        "../envs/samtools.yaml"
    threads: 2
    shell:
        """
        samtools view -h {params.prim_align} {params.minqual} {params.blacklist} \
        {input.bam} {params.atacshift} | 
        samtools view -b > {output} 2> {log}
        """


def get_sambamba_sort_bam(wildcards):
    if sieve_bam(config):
        return rules.sieve_bam.output
    return rules.samtools_presort.output


rule sambamba_sort:
    """
    Sort the result of alignment or sieving with the sambamba sorter.
    """
    input:
        get_sambamba_sort_bam,
    output:
        temp(expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.sambamba-{{sorting}}{{sieve}}.bam", **config)),
    wildcard_constraints:
        sieve="|-sievsort",
    log:
        expand("{log_dir}/sambamba_sort/{{assembly}}-{{sample}}-sambamba_{{sorting}}{{sieve}}.log", **config),
    message: explain_rule("sambamba_sort")
    benchmark:
        expand("{benchmark_dir}/sambamba_sort/{{assembly}}-{{sample}}-{{sorting}}{{sieve}}.benchmark.txt", **config)[0]
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


rule samtools_sort:
    """
    Sort the result of sieving with the samtools sorter.
    """
    input:
        rules.sieve_bam.output,
    output:
        temp(expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-{{sorting}}-sievsort.bam", **config)),
    log:
        expand("{log_dir}/samtools_sort/{{assembly}}-{{sample}}-samtools_{{sorting}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/samtools_sort/{{assembly}}-{{sample}}-{{sorting}}.benchmark.txt", **config)[0]
    params:
        sort_order=lambda wildcards: "-n" if wildcards.sorting == "queryname" else "",
        out_dir=f"{config['result_dir']}/{config['aligner']}",
        memory=lambda wildcards, input, output, threads: f"-m {int(1000 * round(config['bam_sort_mem']/threads, 3))}M",
    threads: 2
    resources:
        mem_gb=config["bam_sort_mem"],
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        # we set this trap to remove temp files when prematurely ending the rule
        trap "rm -f {params.out_dir}/{wildcards.assembly}-{wildcards.sample}.tmp*bam" INT;
        rm -f {params.out_dir}/{wildcards.assembly}-{wildcards.sample}.tmp*bam 2> {log}

        samtools sort {params.sort_order} -@ {threads} {params.memory} {input} -o {output} \
        -T {params.out_dir}/{wildcards.assembly}-{wildcards.sample}.tmp 2> {log}
        """


def get_bam_mark_duplicates(wildcards):
    if sieve_bam(config):
        # when alignmentsieving but not shifting we do not have to re-sort samtools-coordinate
        if wildcards.sorter == "samtools" and wildcards.sorting == "coordinate" and not config.get("tn5_shift", False):
            return rules.sieve_bam.output
        else:
            return expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}-sievsort.bam", **config)

    # if we don't want to do anything get the untreated bam
    if wildcards.sorter == "samtools":
        return rules.samtools_presort.output
    return expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.bam", **config)


rule mark_duplicates:
    """
    Mark all duplicate reads in a bam file with picard MarkDuplicates
    """
    input:
        get_bam_mark_duplicates,
    output:
        bam=expand("{final_bam_dir}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.bam", **config),
        metrics=expand("{qc_dir}/markdup/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.metrics.txt", **config),
    log:
        expand("{log_dir}/mark_duplicates/{{assembly}}-{{sample}}-{{sorter}}-{{sorting}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/mark_duplicates/{{assembly}}-{{sample}}-{{sorter}}-{{sorting}}.benchmark.txt", **config)[0]
    message: explain_rule("mark_duplicates")
    params:
        config["markduplicates"],
    resources:
        mem_gb=5,
    conda:
        "../envs/picard.yaml"
    shell:
        """
        picard MarkDuplicates {params} INPUT={input} \
        OUTPUT={output.bam} METRICS_FILE={output.metrics} > {log} 2>&1
        """


rule samtools_index:
    """
    Create an index of a bam file which can be used for e.g. visualization.
    """
    input:
        "{filepath}.bam",
    output:
        temp("{filepath}.bam.bai") if config.get("cram_no_bam", False) else "{filepath}.bam.bai",
    params:
        config["samtools_index"],
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools index {params} {input} {output}
        """


rule bam2cram:
    """
    Convert a bam file to the more compressed cram format.
    """
    input:
        bam=rules.mark_duplicates.output.bam,
        assembly=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config),
    output:
        expand("{final_bam_dir}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.cram", **config),
    message: explain_rule("bam2cram")
    log:
        expand("{log_dir}/bam2cram/{{assembly}}-{{sample}}-{{sorter}}-{{sorting}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/bam2cram/{{assembly}}-{{sample}}-{{sorter}}-{{sorting}}.benchmark.txt", **config)[0]
    params:
        threads=lambda wildcards, input, output, threads: threads - 1,
    threads: 4
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -T {input.assembly} -C {input.bam} -@ {threads} > {output} 2> {log}
        """


rule samtools_index_cram:
    """
    Generate the index for a cram file.
    """
    input:
        expand("{final_bam_dir}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.cram", **config),
    output:
        expand("{final_bam_dir}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.cram.crai", **config),
    log:
        expand("{log_dir}/samtools_index_cram/{{assembly}}-{{sample}}-{{sorter}}-{{sorting}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/samtools_index_cram/{{assembly}}-{{sample}}-{{sorter}}-{{sorting}}.benchmark.txt", **config)[0]
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools index {input} {output} > {log} 2>&1
        """
