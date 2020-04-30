import os
import re
import gzip


def get_reads(wildcards):
    if 'replicate' in samples and config.get('technical_replicates') == 'merge':
        if config['layout'].get(wildcards.sample, False) == "SINGLE":
            return expand("{trimmed_dir}/merged/{{sample}}_trimmed.{fqsuffix}.gz", **config)
        return sorted(expand("{trimmed_dir}/merged/{{sample}}_{fqext}_trimmed.{fqsuffix}.gz", **config))
    else:
        if config['layout'].get(wildcards.sample, False) == "SINGLE":
            return expand("{trimmed_dir}/{{sample}}_trimmed.{fqsuffix}.gz", **config)
        return sorted(expand("{trimmed_dir}/{{sample}}_{fqext}_trimmed.{fqsuffix}.gz", **config))


if config['aligner'] == 'bowtie2':
    rule bowtie2_index:
        """
        Make a genome index for bowtie2.
        """
        input:
            expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config)
        output:
            directory(expand("{genome_dir}/{{assembly}}/index/bowtie2/", **config))
        log:
            expand("{log_dir}/bowtie2_index/{{assembly}}.log", **config)
        benchmark:
            expand("{benchmark_dir}/bowtie2_index/{{assembly}}.benchmark.txt", **config)[0]
        priority: 1
        threads: 4
        conda:
            "../envs/bowtie2.yaml"
        params: config['index']
        shell:
            "bowtie2-build {params} --threads {threads} {input} {output}/{wildcards.assembly} > {log} 2>&1"


    rule bowtie2_align:
        """
        Align reads against a genome (index) with bowtie2, and pipe the output to the required sorter(s).
        """
        input:
            reads=get_reads,
            index=expand("{genome_dir}/{{assembly}}/index/bowtie2/", **config)
        output:
            pipe(expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-coordinate.pipe", **config)[0])
        log:
            expand("{log_dir}/bowtie2_align/{{assembly}}-{{sample}}.log", **config)
        group: 'alignment'
        benchmark:
            expand("{benchmark_dir}/bowtie2_align/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
        params:
            input=lambda wildcards, input: f'-U {input.reads}' if config['layout'][wildcards.sample] == 'SINGLE' else \
                                           f'-1 {input.reads[0]} -2 {input.reads[1]}',
            params=config['align']
        threads: 20
        conda:
            "../envs/bowtie2.yaml"
        shell:
            """
            bowtie2 {params.params} --threads {threads} -x {input.index}{wildcards.assembly} {params.input} 2> {log} | tee {output} 1> /dev/null 2>> {log}
            """


elif config['aligner'] == 'bwa':
    rule bwa_index:
        """
        Make a genome index for bwa.
        """
        input:
            expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config)
        output:
            directory(expand("{genome_dir}/{{assembly}}/index/bwa/", **config)),
        log:
            expand("{log_dir}/bwa_index/{{assembly}}.log", **config)
        benchmark:
            expand("{benchmark_dir}/bwa_index/{{assembly}}.benchmark.txt", **config)[0]
        params:
            prefix="{genome_dir}/{{assembly}}/index/bwa/{{assembly}}".format(**config),
            params=config['index']
        priority: 1
        resources:
            mem_gb=5
        conda:
            "../envs/bwa.yaml"
        shell:
            "bwa index -p {params.prefix} {params.params} {input} > {log} 2>&1"


    rule bwa_mem:
        """
        Align reads against a genome (index) with bwa, and pipe the output to the required sorter(s).
        """
        input:
            reads=get_reads,
            index=expand("{genome_dir}/{{assembly}}/index/bwa/", **config)
        output:
            pipe(expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-coordinate.pipe", **config)[0])
        log:
            expand("{log_dir}/bwa_mem/{{assembly}}-{{sample}}.log", **config)
        group: 'alignment'
        benchmark:
            expand("{benchmark_dir}/bwa_mem/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
        params:
            index_dir=expand("{genome_dir}/{{assembly}}/index/bwa/{{assembly}}", **config),
            params=config['align']
        resources:
            mem_gb=23
        threads: 20
        conda:
            "../envs/bwa.yaml"
        shell:
            """
            bwa mem {params.params} -t {threads} {params.index_dir} {input.reads} 2> {log} | tee {output} 1> /dev/null 2>> {log}
            """


elif config['aligner'] == 'hisat2':
    rule hisat2_index:
        """
        Make a genome index for hisat2.
        """
        input:
            expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config)
        output:
            directory(expand("{genome_dir}/{{assembly}}/index/hisat2/", **config))
        log:
            expand("{log_dir}/hisat2_index/{{assembly}}.log", **config)
        benchmark:
            expand("{benchmark_dir}/hisat2_index/{{assembly}}.benchmark.txt", **config)[0]
        priority: 1
        threads: 4
        conda:
            "../envs/hisat2.yaml"
        params: config['index']
        shell:
            "hisat2-build {params} -p {threads} {input} {output}/{wildcards.assembly} > {log} 2>&1"


    rule hisat2_align:
        """
        Align reads against a genome (index) with hisat2, and pipe the output to the required sorter(s).
        """
        input:
            reads=get_reads,
            index=expand("{genome_dir}/{{assembly}}/index/hisat2/", **config)
        output:
            pipe(expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-coordinate.pipe", **config)[0])
        log:
            expand("{log_dir}/hisat2_align/{{assembly}}-{{sample}}.log", **config)
        group: 'alignment'
        benchmark:
            expand("{benchmark_dir}/hisat2_align/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
        params:
            input=lambda wildcards, input: f'-U {input.reads}' if config['layout'][wildcards.sample] == 'SINGLE' else \
                                           f'-1 {input.reads[0]} -2 {input.reads[1]}',
            params=config['align']
        threads: 20
        conda:
            "../envs/hisat2.yaml"
        shell:
            """
            hisat2 {params.params} --threads {threads} -x {input.index}{wildcards.assembly} {params.input} 2> {log} | tee {output} 1> /dev/null 2>> {log}
            """


elif config['aligner'] == 'star' or config.get('quantifier', '') == 'star':
    rule star_index:
        """
        Make a genome index for STAR.
        
        Troubleshooting:
        1) sufficient disk space?
        2) increase the RAM available (--limitGenomeGenerateRAM)
        3) reduce the number of threads (snakemake -j 5)
        4) reduce accuracy (--genomeSAsparseD 2)
                
        For example, in your config.yaml, set aligner/quantifier:
        aligner:
            star:
                index: --limitGenomeGenerateRAM 60000000000 --genomeSAsparseD 1
        """
        input:
            genome = expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config),
            sizefile= expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config),
            gtf = expand("{genome_dir}/{{assembly}}/{{assembly}}.gtf", **config)
        output:
            directory(expand("{genome_dir}/{{assembly}}/index/{aligner}", **config))
        log:
            expand("{log_dir}/{aligner}_index/{{assembly}}.log", **config)
        benchmark:
            expand("{benchmark_dir}/{aligner}_index/{{assembly}}.benchmark.txt", **config)[0]
        params:
            config['index']
        priority: 1
        threads: 10
        resources:
            mem_gb=37
        conda:
            "../envs/star.yaml"
        shell:
            """
            function log2 {{
                    local x=0
                    for (( y=$1-1 ; $y > 0; y >>= 1 )) ; do
                        let x=$x+1
                    done
                    echo $x
            }}
            
            # set genome dependent variables
            NBits=""
            NBases=""
            GenomeLength=$(awk -F"\t" '{{x+=$2}}END{{printf "%i", x}}' {input.sizefile})
            NumberOfReferences=$(awk 'END{{print NR}}' {input.sizefile})
            if [ $NumberOfReferences -gt 5000 ]; then
                # for large genomes, --genomeChrBinNbits should be scaled to min(18,log2[max(GenomeLength/NumberOfReferences,ReadLength)])
                # ReadLength is skipped here, as it is unknown
                LpR=$(log2 $((GenomeLength / NumberOfReferences)))
                NBits="--genomeChrBinNbits $(($LpR<18 ? $LpR : 18))"
                printf "NBits: $NBits\n\n" >> {log} 2>&1
            fi
            
            if [ $GenomeLength -lt 268435456 ]; then
                # for small genomes, --genomeSAindexNbases must be scaled down to min(14, log2(GenomeLength)/2-1)
                logG=$(( $(log2 $GenomeLength) / 2 - 1 ))
                NBases="--genomeSAindexNbases $(( $logG<14 ? $logG : 14 ))"
                printf "NBases: $NBases\n\n" >> {log} 2>&1
            fi
            
            mkdir -p {output}
            
            STAR --runMode genomeGenerate --genomeFastaFiles {input.genome} --sjdbGTFfile {input.gtf} \
            --genomeDir {output} --outFileNamePrefix {output}/ \
            --runThreadN {threads} $NBits $NBases {params} >> {log} 2>&1
            """


    rule star_align:
        """
        Align reads against a genome (index) with STAR, and pipe the output to the required sorter(s).
        """
        input:
            reads=get_reads,
            index=expand("{genome_dir}/{{assembly}}/index/{aligner}", **config)
        output:
            dir =directory(expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}", **config)),
            pipe=pipe(expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-coordinate.pipe", **config)[0])
        log:
            directory(expand("{log_dir}/{aligner}_align/{{assembly}}-{{sample}}", **config))
        benchmark:
            expand("{benchmark_dir}/{aligner}_align/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
        params:
            input=lambda wildcards, input: f' {input.reads}' if config['layout'][wildcards.sample] == 'SINGLE' else \
                                           f' {input.reads[0]} {input.reads[1]}',
            params=config['align']
        threads: 8
        resources:
            mem_gb=30
        conda:
            "../envs/star.yaml"
        shell:
            """
            trap "find {log} -type f ! -name Log* -exec rm {{}} \;" EXIT
            mkdir -p {log}
            mkdir -p {output.dir}                

            STAR --genomeDir {input.index} --readFilesIn {params.input} --readFilesCommand gunzip -c \
            --quantMode GeneCounts --outSAMtype BAM Unsorted --outStd BAM_Unsorted \
            --outFileNamePrefix {log}/ --outTmpDir {output.dir}/STARtmp \
            --runThreadN {threads} {params.params} > {output.pipe} 2> {log}/Log.stderr.out

            # move all non-log files to output directory (this way the log files are kept on error)
            find {log} -type f ! -name Log* -exec mv {{}} {output.dir} \;
            """


rule samtools_presort:
    """
    Sort the result of alignment with the samtools sorter.
    """
    input:
        expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-coordinate.pipe", **config)
    output:
        temp(expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-coordinate-unsieved.bam", **config))
    log:
        expand("{log_dir}/samtools_presort/{{assembly}}-{{sample}}.log", **config)
    group: 'alignment'
    benchmark:
        expand("{benchmark_dir}/samtools_presort/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
    params:
        threads=lambda wildcards, input, output, threads: max([1, threads - 1]),
        sort_order="-n" if not use_alignmentsieve(config) and config.get("bam_sort_order") == "queryname" else "",
        out_dir=f"{config['result_dir']}/{config['aligner']}"
    threads: 3
    resources:
        mem_mb=2500
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        trap "rm -f {params.out_dir}/{wildcards.assembly}-{wildcards.sample}.tmp*bam" INT;
        samtools sort {params.sort_order} -@ {params.threads} {input} -o {output} -T {params.out_dir}/{wildcards.assembly}-{wildcards.sample}.tmp 2> {log}
        """


def get_blacklist_files(wildcards):
    files = {}
    # ideally get genome is a checkpoint, however there are quite some Snakemake
    # bugs related to this. So for now we solve it like this
    # TODO: switch back to checkpoints
    if config.get('remove_blacklist') and wildcards.assembly.lower() in \
            ["ce10", "dm3", "hg38", "hg19", "mm9", "mm10"]:
        blacklist = f"{config['genome_dir']}/{wildcards.assembly}/{wildcards.assembly}.fa"
        files['blacklist'] = blacklist

    if config.get('remove_mito'):
        sizes = f"{config['genome_dir']}/{wildcards.assembly}/{wildcards.assembly}.fa.sizes"
        files['sizes'] = sizes

    return files


rule setup_blacklist:
    input:
        unpack(get_blacklist_files)
    output:
        temp(expand("{genome_dir}/{{assembly}}/{{assembly}}.customblacklist.bed", **config))
    log:
        expand("{log_dir}/setup_blacklist/{{assembly}}.log", **config)
    run:
        newblacklist = ""
        if config.get('remove_blacklist') and wildcards.assembly.lower() in \
                ["ce10", "dm3", "hg38", "hg19", "mm9", "mm10"]:
            blacklist = f"{config['genome_dir']}/{wildcards.assembly}/{wildcards.assembly}.blacklist.bed.gz"
            with gzip.GzipFile(blacklist) as file:
                newblacklist += file.read().decode('utf8')

        if any('.fa.sizes' in inputfile for inputfile in input):
            with open(input.sizes, 'r') as file:
                sizesfile = file.read().strip()
                for match in re.findall("chrM.*|chrm.*|MT.*", sizesfile):
                    chrm, size = match.split('\t')
                    newblacklist += f"{chrm}\t0\t{size}\n"

        with open(output[0], 'w') as f:
            f.write(newblacklist)


rule alignmentsieve:
    input:
        bam=rules.samtools_presort.output,
        bai=expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-coordinate-unsieved.bam.bai", **config),
        blacklist=rules.setup_blacklist.output
    output:
        temp(expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-coordinate-sieved.bam", **config))
    log:
        expand("{log_dir}/alignmentsieve/{{assembly}}-{{sample}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/alignmentsieve/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
    params:
        minqual=f"--minMappingQuality {config['min_mapping_quality']}",
        atacshift="--ATACshift" if config['tn5_shift'] else '',
        blacklist=lambda wildcards, input: "" if os.path.exists(input.blacklist[0]) and os.stat(input.blacklist[0]).st_size == 0 else \
                                          f"--blackListFileName {input.blacklist}"
    conda:
        "../envs/deeptools.yaml"
    threads: 4
    resources:
        deeptools_limit=lambda wildcards, threads: threads
    shell:
        "alignmentSieve -b {input.bam} -o {output} {params.minqual} {params.atacshift} {params.blacklist} -p {threads} > {log} 2>&1"


def get_sambamba_sort_bam(wildcards):
    if use_alignmentsieve(config):
        return rules.alignmentsieve.output
    return rules.samtools_presort.output


rule sambamba_sort:
    """
    Sort the result of alignment or sieving with the sambamba sorter.
    """
    input:
        get_sambamba_sort_bam
    output:
        temp(expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.sambamba-{{sorting}}{{sieve}}.bam", **config))
    wildcard_constraints:
        sieve="|-sievsort"
    log:
        expand("{log_dir}/sambamba_sort/{{assembly}}-{{sample}}-sambamba_{{sorting}}{{sieve}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/sambamba_sort/{{assembly}}-{{sample}}-{{sorting}}{{sieve}}.benchmark.txt", **config)[0]
    params:
        lambda wildcards: "-n" if 'queryname' in wildcards.sorting else '',
    threads: 2
    conda:
        "../envs/sambamba.yaml"
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
        rules.alignmentsieve.output
    output:
        temp(expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-{{sorting}}-sievsort.bam", **config))
    log:
        expand("{log_dir}/samtools_sort/{{assembly}}-{{sample}}-samtools_{{sorting}}.log", **config)
    group: 'alignment'
    benchmark:
        expand("{benchmark_dir}/samtools_sort/{{assembly}}-{{sample}}-{{sorting}}.benchmark.txt", **config)[0]
    params:
        order=lambda wildcards: "-n" if wildcards.sorting == 'queryname' else '',
        threads=lambda wildcards, input, output, threads: max([1, threads - 1]),
        out_dir=f"{config['result_dir']}/{config['aligner']}"
    threads: 2
    resources:
        mem_mb=2500
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        trap "rm -f {params.out_dir}/{wildcards.assembly}-{wildcards.sample}.tmp*bam" INT;
        samtools sort {params.sort_order} -@ {params.threads} {input} -o {output} -T {params.out_dir}/{wildcards.assembly}-{wildcards.sample}.tmp 2> {log}
        """


def get_bam_mark_duplicates(wildcards):
    if use_alignmentsieve(config):
        # when alignmentsieving but not shifting we do not have to re-sort samtools-coordinate
        if wildcards.sorter == 'samtools' and wildcards.sorting == 'coordinate' and \
                not config.get('tn5_shift', False):
            return rules.alignmentsieve.output
        else:
            return expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}-sievsort.bam", **config)

    # if we don't want to do anything get the untreated bam
    if wildcards.sorter == 'samtools':
        return rules.samtools_presort.output
    return expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.bam", **config)


rule mark_duplicates:
    """
    Mark (but keep) all duplicate reads in a bam file with picard MarkDuplicates
    """
    input:
        get_bam_mark_duplicates
    output:
        bam=    expand("{dedup_dir}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.bam", **config),
        metrics=expand("{qc_dir}/dedup/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.metrics.txt", **config)
    log:
        expand("{log_dir}/mark_duplicates/{{assembly}}-{{sample}}-{{sorter}}-{{sorting}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/mark_duplicates/{{assembly}}-{{sample}}-{{sorter}}-{{sorting}}.benchmark.txt", **config)[0]
    params:
        config['markduplicates']
    resources:
        mem_gb=10
    conda:
        "../envs/picard.yaml"
    shell:
        "picard MarkDuplicates {params} INPUT={input} "
        "OUTPUT={output.bam} METRICS_FILE={output.metrics} > {log} 2>&1"


rule samtools_index:
    """
    Create an index of a bam file which can be used for e.g. visualization.
    """
    input:
        "{filepath}.bam"
    output:
        temp("{filepath}.bam.bai") if config.get("cram_no_bam", False) else "{filepath}.bam.bai"
    params:
        config['samtools_index']
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools index {params} {input} {output}
        """

        
rule bam2cram:
    """
    Convert bam to the more compressed cram format
    """
    input:
         bam=rules.mark_duplicates.output.bam,
         assembly=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config)
    output:
        expand("{dedup_dir}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.cram", **config),
    log:
        expand("{log_dir}/bam2cram/{{assembly}}-{{sample}}-{{sorter}}-{{sorting}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/bam2cram/{{assembly}}-{{sample}}-{{sorter}}-{{sorting}}.benchmark.txt", **config)[0]
    params:
        threads=lambda wildcards, input, output, threads: threads - 1
    threads: 4
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -@ {threads} -T {input.assembly} -C {input.bam} > {output} 2> {log}"


rule samtools_index_cram:
    input:
        expand("{dedup_dir}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.cram", **config),
    output:
        expand("{dedup_dir}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.cram.crai", **config),
    log:
        expand("{log_dir}/samtools_index_cram/{{assembly}}-{{sample}}-{{sorter}}-{{sorting}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/samtools_index_cram/{{assembly}}-{{sample}}-{{sorter}}-{{sorting}}.benchmark.txt", **config)[0]
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input} {output} > {log} 2>&1"
