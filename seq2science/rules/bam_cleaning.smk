"""
All rules/logic related to filtering (sieving) after alignment to a genome should be here.
"""

from seq2science.util import sieve_bam


localrules:
    setup_blacklist,
    complement_blacklist,


# the blacklist has different output depending on what it blacklists
blacklisted_filename = []
if config.get("remove_blacklist"):
    blacklisted_filename.append("encode")
if config.get("remove_mito"):
    blacklisted_filename.append("mito")
blacklisted_filename = "_".join(blacklisted_filename)


rule setup_blacklist:
    """
    Combine the encode blacklist with mitochrondrial dna depending on config.
    """
    input:
        blacklist=expand("{genome_dir}/{{assembly}}/{{assembly}}.blacklist.bed", **config),
        sizes=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config),
    output:
        expand(
            "{genome_dir}/{{assembly}}/{{assembly}}.seq2scienceblacklist_{types}.bed",
            **{**config, **{"types": blacklisted_filename}}
        ),
    params:
        config.get("remove_blacklist"),  # helps resolve changed params
        config.get("remove_mito"),  # helps resolve changed params
    run:
        newblacklist = ""
        if config.get("remove_blacklist"):
            with open(input["blacklist"][0]) as file:
                newblacklist += file.read()

        if config.get("remove_mito"):
            with open(input["sizes"][0]) as file:
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
        expand(
            "{genome_dir}/{{assembly}}/{{assembly}}.seq2scienceblacklist_complement_{types}.bed",
            **{**config, **{"types": blacklisted_filename}}
        ),
    params:
        config.get("remove_blacklist"),  # helps resolve changed params
        config.get("remove_mito"),  # helps resolve changed params
    log:
        expand("{log_dir}/complement_blacklist/{{assembly}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/complement_blacklist/{{assembly}}.benchmark.txt", **config)[0]
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        sortBed -faidx {input.sizes} -i {input.blacklist} 2>> {log} |
        complementBed -i stdin -g {input.sizes} > {output} 2>> {log}
        """


rule mark_duplicates:
    """
    Mark all duplicate reads in a bam file with picard MarkDuplicates
    """
    input:
        rules.samtools_presort.output,
    output:
        bam=temp(expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-coordinate-dupmarked.bam", **config)),
        metrics=expand("{qc_dir}/markdup/{{assembly}}-{{sample}}.samtools-coordinate.metrics.txt", **config),
    log:
        expand("{log_dir}/mark_duplicates/{{assembly}}-{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/mark_duplicates/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
    message:
        explain_rule("mark_duplicates")
    params:
        config["markduplicates"],
    resources:
        mem_gb=8,
    conda:
        "../envs/picard.yaml"
    wildcard_constraints:
        sample=f"""({any_given("sample", "technical_replicates", "control")})(_allsizes)?""",
    shell:
        """
        # use the TMPDIR if set, and not given in the config
        if [[ ${{TMPDIR:=F}} == "F" ]] || [[ "{params}" == *TMP_DIR* ]]
        then
            tmpdir=""
        else
            tmpdir=TMP_DIR=$TMPDIR
        fi
        picard MarkDuplicates $tmpdir {params} INPUT={input} \
        OUTPUT={output.bam} METRICS_FILE={output.metrics} > {log} 2>&1
        """


# if doing tn5 shift, we need to re-sort afterwards!
if config.get("tn5_shift"):
    shiftsieve = "-shifted"
    sieve_bam_output = {"final": temp(f"{config['final_bam_dir']}/{{assembly}}-{{sample}}.samtools-coordinate{shiftsieve}.bam")}
else:
    shiftsieve = ""
    sieve_bam_output = {"final": f"{config['final_bam_dir']}/{{assembly}}-{{sample}}.samtools-coordinate.bam"}

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
        sizes=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config),
    output:
        **sieve_bam_output,
    log:
        expand("{log_dir}/sieve_bam/{{assembly}}-{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/sieve_bam/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
    message:
        explain_rule("sieve_bam")
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
            if sampledict[wildcards.sample]["layout"] == "PAIRED" and config["filter_on_size"]
            else ""
        ),
        sizesieve_touch=(
            lambda wildcards, input, output: f"touch {output.allsizes}"
            if sampledict[wildcards.sample]["layout"] == "SINGLE" and config["filter_on_size"]
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


if not sieve_bam(config):

    rule cp_unsieved2sieved:
        """
        Copy a bam file if no sieving is necessary. 
        """
        input:
            rules.mark_duplicates.output.bam,
        output:
            expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config),
        shell:
            """
            cp {input} {output}
            """

    ruleorder: cp_unsieved2sieved > sieve_bam > samtools_sort


def get_sambamba_sort_bam(wildcards):
    if sieve_bam(config):
        return rules.sieve_bam.output
    return rules.samtools_presort.output


rule sambamba_sort:
    """
    Sort the result of alignment or sieving with the sambamba sorter.
    """
    input:
        expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config),
    output:
        temp(expand("{final_bam_dir}/{{assembly}}-{{sample}}.sambamba-{{sorting}}.bam", **config)),
    wildcard_constraints:
        sieve="|-sievsort",
    log:
        expand("{log_dir}/sambamba_sort/{{assembly}}-{{sample}}-sambamba_{{sorting}}.log", **config),
    message:
        explain_rule("sambamba_sort")
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


rule samtools_sort:
    """
    Sort the result of shiftsieving with the samtools sorter.
    """
    input:
        rules.sieve_bam.output.final,
    output:
        expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-{{sorting}}.bam", **config),
    log:
        expand("{log_dir}/samtools_sort/{{assembly}}-{{sample}}-samtools_{{sorting}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/samtools_sort/{{assembly}}-{{sample}}-{{sorting}}.benchmark.txt", **config)[0]
    params:
        sort_order=lambda wildcards: "-n" if wildcards.sorting == "queryname" else "",
        out_dir=f"{config['result_dir']}/{config['aligner']}",
        memory=config['bam_sort_mem'],
    wildcard_constraints:
        sample=f"""({any_given("sample", "technical_replicates", "control")})""",
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

        # RAM per thread in MB
        memory=$((1000*{params.memory}/{threads}))M

        samtools sort {params.sort_order} -@ {threads} -m $memory {input} -o {output} \
        -T {params.out_dir}/{wildcards.assembly}-{wildcards.sample}.tmp 2> {log}
        """


if config["filter_on_size"]:
    rule samtools_sort_allsizes:
        """
        Sort the result of shiftsieving with the samtools sorter. This rule is identical to the
        samtools_sort rule except that the output is temporary.
        """
        input:
            rules.sieve_bam.output.final,
        output:
            temp(expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-{{sorting}}.bam", **config)),
        log:
            expand("{log_dir}/samtools_sort/{{assembly}}-{{sample}}-samtools_{{sorting}}.log", **config),
        benchmark:
            expand("{benchmark_dir}/samtools_sort/{{assembly}}-{{sample}}-{{sorting}}.benchmark.txt", **config)[0]
        params:
            sort_order=lambda wildcards: "-n" if wildcards.sorting == "queryname" else "",
            out_dir=f"{config['result_dir']}/{config['aligner']}",
            memory=config['bam_sort_mem'],
        wildcard_constraints:
            sample=f"""({any_given("sample", "technical_replicates", "control")})_allsizes""",
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
            # RAM per thread in MB
            memory=$((1000*{params.memory}/{threads}))M
            samtools sort {params.sort_order} -@ {threads} -m $memory {input} -o {output} \
            -T {params.out_dir}/{wildcards.assembly}-{wildcards.sample}.tmp 2> {log}
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


if config["filter_on_size"]:

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


rule bam2cram:
    """
    Convert a bam file to the more compressed cram format.
    """
    input:
        bam=rules.mark_duplicates.output.bam,
        assembly=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config),
    output:
        expand("{final_bam_dir}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.cram", **config),
    message:
        explain_rule("bam2cram")
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
        expand(
            "{benchmark_dir}/samtools_index_cram/{{assembly}}-{{sample}}-{{sorter}}-{{sorting}}.benchmark.txt", **config
        )[0]
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools index {input} {output} > {log} 2>&1
        """
