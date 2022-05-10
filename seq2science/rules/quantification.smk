"""
all rules/logic related to counting/quantification of genes should be here.
"""

import os
import os.path

from seq2science.util import get_bustools_rid

if config["quantifier"] == "salmon":

    rule get_transcripts:
        """
        Generate transcripts.fasta using gffread.

        Requires genome.fa and annotation.gtf (with matching chromosome/scaffold names)
        """
        input:
            fa=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config),
            gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config),
        output:
            expand("{genome_dir}/{{assembly}}/{{assembly}}.transcripts.fa", **config),
        log:
            expand("{log_dir}/quantification/{{assembly}}.transcripts.log", **config),
        conda: "../envs/salmon.yaml"
        priority: 1
        shell:
            """
            gffread -w {output} -g {input.fa} {input.gtf} > {log} 2>&1
            """

    rule partial_decoy_transcripts:
        """
        Compute a set of decoy sequences by mapping the annotated transcripts against 
        a hard-masked version of the organism’s genome.
        """
        input:
            genome=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config),
            gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config),
            transcripts=rules.get_transcripts.output,
        output:
            gentrome=temp(expand("{genome_dir}/{{assembly}}/gentrome.fa", **config)),
            decoys=temp(expand("{genome_dir}/{{assembly}}/decoys.txt", **config)),
        params:
            script=f"{config['rule_dir']}/../scripts/generateDecoyTranscriptome.sh",
        log:
            expand("{log_dir}/quantification/{{assembly}}.partial_decoy_transcripts.log",**config),
        message: EXPLAIN["partially_decoy_aware"]
        conda: "../envs/decoy.yaml"
        threads: 40
        resources: mem_gb=64,
        priority: 1
        shell:
            ("cpulimit --include-children -l {threads}00 --" if config.get("cpulimit", True) else "")+
            """\
            sh {params.script} -j {threads} -g {input.genome} -a {input.gtf} -t {input.transcripts} -o $(dirname {output[0]}) > {log} 2>&1
            """

    rule full_decoy_transcripts:
        """
        Use the entire genome as decoy sequence.
        """
        input:
            genome=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa",**config),
            sizes=rules.get_genome_support_files.output.sizes,
            transcripts=rules.get_transcripts.output,
        output:
            gentrome=temp(expand("{genome_dir}/{{assembly}}/full_gentrome.fa",**config)),
            decoys=temp(expand("{genome_dir}/{{assembly}}/full_decoys.txt",**config)),
        message: EXPLAIN["fully_decoy_aware"]
        shell:
            """
            cut -f 1 {input.sizes} > {output.decoys}
            cat {input.transcripts} {input.genome} > {output.gentrome}
            """


    def salmon_index_input(wildcards):
        level = config.get("quantifier_decoys", "")
        if level == "full":
            return rules.full_decoy_transcripts.output
        if level == "partial":
            return rules.partial_decoy_transcripts.output
        return rules.get_transcripts.output


    salmon_index_output = "salmon"
    if config.get("quantifier_decoys", "") == "full":
        salmon_index_output += "_fully_decoy_aware"
    elif config.get("quantifier_decoys", "") == "partial":
        salmon_index_output += "_partially_decoy_aware"


    rule salmon_index:
        """
        Generate a transcriptome index for Salmon.
        """
        input:
            salmon_index_input
        output:
            directory(expand(f"{{genome_dir}}/{{{{assembly}}}}/index/{salmon_index_output}", **config)),
        log:
            expand(f"{{log_dir}}/quantification/{{{{assembly}}}}_{salmon_index_output}_index.log",**config),
        params:
            decoys=(
                lambda wildcards, input: [""]
                if config.get("quantifier_decoys", "none") == "none"
                else ["-d", input[1]]
            ),
            params=config["quantifier_index"],
        conda: "../envs/salmon.yaml"
        priority: 1
        threads: 10
        resources: mem_gb=9,  # fully decoy aware
        shell:
            """
            mkdir -p {output}

            salmon index -t {input[0]} {params.decoys} -i {output} {params.params} \
            --threads {threads} > {log} 2>&1
            """

    rule salmon_quant:
        """
        Align reads against a transcriptome (index) with Salmon (mapping-based mode) and output a quantification file per sample.
        """
        input:
            reads=get_reads,
            index=rules.salmon_index.output,
        output:
            expand("{result_dir}/{quantifier}/{{assembly}}-{{sample}}/quant.sf", **config),
        log:
            expand("{log_dir}/quantification/{{assembly}}_{quantifier}-{{sample}}.log",**config),
        params:
            input=(
                lambda wildcards, input: ["-r", input.reads]
                if SAMPLEDICT[wildcards.sample]["layout"] == "SINGLE"
                else ["-1", input.reads[0], "-2", input.reads[1]]
            ),
            params=config["quantifier_flags"],
            reps=lambda wildcards, input: input,  # help resolve changes in input files
        message: EXPLAIN["salmon_quant"]
        conda: "../envs/salmon.yaml"
        threads: 12
        resources: mem_gb=8,
        shell:
            """
            salmon quant -i {input.index} -l A {params.input} {params.params} -o $(dirname {output}) \
            --threads {threads} > {log} 2>&1
            """

elif  "scrna_seq" == WORKFLOW:

    def get_fastq_pair_reads(wildcards):
        """
        Extracts the correct combination of R1/R2 (trimmed and barcodes) for fastq_pair
        based on Kallisto bustools settings.
        """
        reads = dict()
        assert (
            SAMPLEDICT[wildcards.sample]["layout"] == "PAIRED"
        ), "Seq2science does not support scRNA-seq samples that are single-ended"

        if config["quantifier"] == "kallistobus":
            read_id = get_bustools_rid(config.get("count"))
            # Determine mate for trimming
            if read_id == 0:
                reads["r1"] = expand("{trimmed_dir}/{{sample}}_{fqext1}_trimmed.{fqsuffix}.gz", **config)
                reads["r2"] = expand("{fastq_dir}/{{sample}}_{fqext2}.{fqsuffix}.gz", **config)
            elif read_id == 1:
                reads["r1"] = expand("{fastq_dir}/{{sample}}_{fqext1}.{fqsuffix}.gz", **config)
                reads["r2"] = expand("{trimmed_dir}/{{sample}}_{fqext2}_trimmed.{fqsuffix}.gz", **config)

        elif config["quantifier"] == "citeseqcount":
            reads["r1"] = expand("{fastq_dir}/{{sample}}_{fqext1}.{fqsuffix}.gz", **config)
            reads["r2"] = expand("{trimmed_dir}/{{sample}}_{fqext2}_trimmed.{fqsuffix}.gz", **config)
        else:
            logger.error(
                f"Something went wrong parsing the read id for fastq_pair. "
                "Please make an issue on github if this is unexpected behaviour!"
            )
            os._exit(0)  # noqa

        return reads


    rule fastq_pair:
        """
        fastq_pair re-writes paired-end fastq files to ensure that each read has a mate and
        discards singleton reads. This step is required after scRNA trimming since we only trim the fastq
        containing reads and not the barcode fastq.
        """
        input:
            unpack(get_fastq_pair_reads),
        output:
            reads=temp(expand("{fastq_clean_dir}/{{sample}}_clean_{fqext}.{fqsuffix}.paired.fq.gz", **config)),
            intermediates1=temp(expand("{fastq_clean_dir}/{{sample}}_clean_{fqext}.{fqsuffix}", **config)),
            intermediates2=temp(expand("{fastq_clean_dir}/{{sample}}_clean_{fqext}.{fqsuffix}.single.fq", **config)),
            intermediates3=temp(expand("{fastq_clean_dir}/{{sample}}_clean_{fqext}.{fqsuffix}.paired.fq", **config)),
        priority: 1
        log:
            expand("{log_dir}/fastq_pair/{{sample}}.log", **config),
        benchmark:
            expand("{benchmark_dir}/fastq_pair/{{sample}}.benchmark.txt", **config)[0]
        conda:
            "../envs/fastq-pair.yaml"
        params:
            options=config.get("fastq-pair", ""),
            tused=lambda wildcards, input: "true" if "-t" in config.get("fastq-pair", "") else "false",
        shell:
            """
            gunzip -c {input.r1} > {output.intermediates1[0]} 2> {log}
            gunzip -c {input.r2} > {output.intermediates1[1]} 2>> {log}
            if [ {params.tused} == true ]
            then
              opts="{params.options}"
            else
              echo "\nsetting parameter t with the number of reads in the fastq\n" >> {log}
              opts="-p -t "$(wc -l {input.r1} | grep -Po '^\d+' | awk '{{print int($1/4)}}')
            fi
            fastq_pair $opts {output.intermediates1} >> {log} 2>&1
            gzip -c {output.intermediates3[0]} > {output.reads[0]} 2> {log}
            gzip -c {output.intermediates3[1]} > {output.reads[1]} 2>> {log}
            """


    def get_final_reads(wildcards):
        if wildcards.sample in MERGED_TREPS:
            return expand(f"{{trimmed_dir}}/{wildcards.sample}_{{fqext}}_trimmed.{{fqsuffix}}.gz", **config)
        return rules.fastq_pair.output.reads


    if config["quantifier"] == "kallistobus":
        if "kite" in config.get("ref", ""):

            ruleorder: kallistobus_ref_kite > get_genome
            ruleorder: kallistobus_ref_kite > kallistobus_ref


        else:

            ruleorder: kallistobus_ref > kallistobus_ref_kite

        rule kallistobus_ref:
            """
            Make a genome index for kallistobus. This index is required for counting.
            """
            input:
                fa=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config),
                gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config),
            output:
                directory(expand("{genome_dir}/{{assembly}}/index/kallistobus", **config)),
            log:
                expand("{log_dir}/kallistobus_index/{{assembly}}.log", **config),
            benchmark:
                expand("{benchmark_dir}/kallistobus_index/{{assembly}}.benchmark.txt", **config)[0]
            priority: 1
            conda:
                "../envs/kallistobus.yaml"
            resources:
                mem_gb=88,
            params:
                basename=lambda wildcards, output: f"{output[0]}/{wildcards.assembly}",
                options=config.get("ref"),
                c1=lambda wildcards, output: f"-c1 {output[0]}/{wildcards.assembly}_cdna_t2c.txt"
                if ("lamanno" or "nucleus") in config.get("ref")
                else "",
                c2=lambda wildcards, output: f"-c2 {output[0]}/{wildcards.assembly}_intron_t2c.txt"
                if ("lamanno" or "nucleus") in config.get("ref")
                else "",
                f2=lambda wildcards, output: f"-f2 {output[0]}/{wildcards.assembly}_intron.fa"
                if ("lamanno" or "nucleus") in config.get("ref")
                else "",
            shell:
                """
                mkdir -p {output}
                kb ref \
                {input.fa} {input.gtf} \
                -i {params.basename}.idx -g {params.basename}_t2g.txt -f1 {params.basename}_cdna.fa \
                {params.f2} {params.c1} {params.c2} \
                {params.options} > {log} 2>&1
                """
                    
                    
        rule kallistobus_ref_kite:
            """
            Make a mismatch index for kallistobus. This index is required to count feature barcodes, such as antibody tags.
            """
            input:
                featurebarcodes=expand("{genome_dir}/{{assembly}}.tsv", **config),
            output:
                directory(expand("{genome_dir}/{{assembly}}/index/kallistobus/kite", **config)),
            log:
                expand("{log_dir}/kallistobus_index_kite/{{assembly}}.log", **config),
            conda:
                "../envs/kallistobus.yaml"
            resources:
                mem_gb=12,
            params:
                basename=lambda wildcards, output: f"{output[0]}/{wildcards.assembly}",
                options=config.get("ref"),
            priority: 1
            shell:
                """
                mkdir -p {output}
                kb ref  \
                {input.featurebarcodes} \
                {params.options} \
                -i {params.basename}.idx -g {params.basename}_t2g.txt -f1 {params.basename}_cdna.fa > {log} 2>&1
                """


        def get_kb_dir(wildcards):
            # Get the correct assembly directory for ADT/Genome based workflows
            if "kite" in config.get("ref", ""):
                return directory(expand("{genome_dir}/{{assembly}}/index/kallistobus/kite", **config))
            else:
                return directory(expand("{genome_dir}/{{assembly}}/index/kallistobus", **config))

        rule kallistobus_count:
            """
            Align reads against a transcriptome (index) with kallistobus and output a quantification file per sample.
            """
            input:
                basedir=get_kb_dir,
                reads=get_final_reads,
                barcodefile=config.get("barcodefile", []),
            output:
                dir=expand(
                    "{result_dir}/{quantifier}/{{assembly}}-{{sample}}/{file}",
                    **{**config, **{"file": ["inspect.json", "run_info.json", "output.bus"]}}
                ),
            log:
                expand("{log_dir}/kallistobus_count/{{assembly}}-{{sample}}.log", **config),
            benchmark:
                expand("{benchmark_dir}/kallistobus_count/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
            priority: 1
            conda:
                "../envs/kallistobus.yaml"
            threads: 8
            message: EXPLAIN["kallistobus-count"]
            resources:
                mem_gb=66,
            params:
                basename=lambda wildcards, input: f"{input.basedir[0]}/{wildcards.assembly}",
                barcode_arg=lambda wildcards, input: ("-w " + input.barcodefile) if input.barcodefile else "",
                options=config.get("count"),
                outdir=lambda wildcards, input, output: os.path.dirname(output[0]),
                c1=lambda wildcards, input: f"-c1 {input[0]}/{wildcards.assembly}_cdna_t2c.txt"
                if ("lamanno" or "nucleus") in config.get("count")
                else "",
                c2=lambda wildcards, input: f"-c2 {input[0]}/{wildcards.assembly}_intron_t2c.txt"
                if ("lamanno" or "nucleus") in config.get("count")
                else "",
            shell:
                """
                kb count \
                -i {params.basename}.idx \
                -t {threads} -g {params.basename}_t2g.txt \
                -o {params.outdir} {params.c1} {params.c2} \
                {params.barcode_arg} {params.options} {input.reads} > {log} 2>&1
                # Validate output
                if grep -q 'ERROR\|bad_alloc' "{log}"; then
                  exit 1
                fi
                """


    if config["quantifier"] == "citeseqcount":

        ruleorder: citeseqcount > get_genome

        rule citeseqcount:
            """
            ADT mapping and quantification with cite-seq-count and output a umi/read matrix per sample.
            """
            input:
                reads=get_final_reads,
                barcodefile=config.get("barcodefile", []),
                tags=expand("{genome_dir}/{{assembly}}.csv", **config),
            output:
                dir=expand(
                    "{result_dir}/{quantifier}/{{assembly}}-{{sample}}/_tag_counts/{file}",
                    **{**config, **{"file": ["run_report.yaml", "unmapped.csv"]}}
                ),
            log:
                expand("{log_dir}/citeseqcount/{{assembly}}-{{sample}}.log", **config),
            benchmark:
                expand("{benchmark_dir}/citeseqcount/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
            priority: 1
            conda:
                "../envs/cite-seq-count.yaml"
            threads: 8
            message: EXPLAIN["citeseqcount"]
            resources:
                mem_gb=66,
            params:
                barcode_arg=lambda wildcards, input: ("-wl " + input.barcodefile) if input.barcodefile else "",
                options=config.get("count"),
                outdir=lambda wildcards, input, output: os.path.dirname(output[0]),
            shell:
                """
                CITE-seq-Count -R1 {input.reads[0]} -R2 {input.reads[1]} \
                -t {input.tags} {params.barcode_arg} -T {threads} \
                {params.options} -o {params.outdir} > {log} 2>&1
                """

elif config["quantifier"] == "htseq":

    rule htseq_count:
        """
        summarize reads to gene level. Outputs a counts table per bam file.
        """
        input:
            bam=expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam",**config),
            gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf",**config),
            report=rules.infer_strandedness.output,
        output:
            expand("{counts_dir}/{{assembly}}-{{sample}}.counts.tsv",**config),
        log:
            expand("{log_dir}/counts_matrix/{{assembly}}-{{sample}}.counts.log",**config),
        params:
            strandedness=lambda wildcards, input: get_strandedness(input.report[0]),
            user_flags=config["htseq_flags"],
        message: EXPLAIN["htseq_count"]
        conda: "../envs/gene_counts.yaml"
        threads: 1
        shell:
            """
            htseq-count {input.bam} {input.gtf} -r pos -s {params.strandedness} {params.user_flags} -n {threads} -c {output} > {log} 2>&1
            """

elif config["quantifier"] == "featurecounts":

    rule featurecounts:
        """
        summarize reads to gene level. Outputs a counts table per bam file.
        """
        input:
            bam=expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam",**config),
            gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf",**config),
            report=rules.infer_strandedness.output,
        output:
            expand("{counts_dir}/{{assembly}}-{{sample}}.counts.tsv",**config),
        log:
            expand("{log_dir}/counts_matrix/{{assembly}}-{{sample}}.counts.log",**config),
        params:
            strandedness=lambda wildcards, input: get_strandedness(input.report[0], fmt="fc"),
            endedness=lambda wildcards: "" if SAMPLEDICT[wildcards.sample]["layout"] == "SINGLE" else "-p",
            user_flags=config["featurecounts_flags"],
        message: EXPLAIN["featurecounts_rna"]
        conda: "../envs/gene_counts.yaml"
        threads: 1
        shell:
            """
            featureCounts -a {input.gtf} {input.bam} {params.endedness} -s {params.strandedness} {params.user_flags} -T {threads} -o {output} > {log} 2>&1
            """

if config.get("dexseq"):

    rule prepare_DEXseq_annotation:
        """
        generate a DEXseq annotation.gff from the annotation.gtf
        """
        input:
            expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf",**config),
        output:
            expand("{genome_dir}/{{assembly}}/{{assembly}}.DEXseq_annotation.gff",**config),
        log:
            expand("{log_dir}/counts_matrix/{{assembly}}.prepare_DEXseq_annotation.log",**config),
        conda: "../envs/dexseq.yaml"
        shell:
            """
            current_conda_env=$(conda env list | grep \* | cut -d "*" -f2-)
            DEXseq_path=${{current_conda_env}}/lib/R/library/DEXSeq/python_scripts

            python ${{DEXseq_path}}/dexseq_prepare_annotation.py {input} {output} > {log} 2>&1
            """

    rule dexseq_count:
        """
        count exon usage
        """
        input:
            bam=expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam",**config),
            bai=expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam.bai",**config),
            gff=expand("{genome_dir}/{{assembly}}/{{assembly}}.DEXseq_annotation.gff",**config),
            report=rules.infer_strandedness.output,
        output:
            expand("{counts_dir}/{{assembly}}-{{sample}}.DEXSeq_counts.tsv",**config),
        log:
            expand("{log_dir}/counts_matrix/{{assembly}}-{{sample}}.DEXseq_counts.log",**config),
        params:
            strandedness=lambda wildcards, input: get_strandedness(input.report[0]),
            endedness=lambda wildcards: "" if SAMPLEDICT[wildcards.sample]["layout"] == "SINGLE" else "-p yes",
        message: EXPLAIN["dexseq"]
        conda: "../envs/dexseq.yaml"
        shell:
            """
            current_conda_env=$(conda env list | grep \* | cut -d "*" -f2-)
            DEXseq_path=${{current_conda_env}}/lib/R/library/DEXSeq/python_scripts

            python ${{DEXseq_path}}/dexseq_count.py -f bam -r pos {params.endedness} -s {params.strandedness} {input.gff} {input.bam} {output} > {log} 2>&1
            """
