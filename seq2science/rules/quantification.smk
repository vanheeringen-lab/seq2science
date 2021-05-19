import os
import glob

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
            expand("{log_dir}/get_genome/{{assembly}}.transcripts.log", **config),
        benchmark:
            expand("{benchmark_dir}/get_genome/{{assembly}}.transcripts.benchmark.txt", **config)[0]
        conda:
            "../envs/salmon.yaml"
        priority: 1
        shell:
            "gffread -w {output} -g {input.fa} {input.gtf} >> {log} 2>&1"


    rule decoy_transcripts:
        """
        Generate decoy_transcripts.txt for Salmon indexing  
    
        script source: https://github.com/COMBINE-lab/SalmonTools
        """
        input:
            genome=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config),
            gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config),
            transcripts=expand("{genome_dir}/{{assembly}}/{{assembly}}.transcripts.fa", **config),
        output:
            gentrome=expand("{genome_dir}/{{assembly}}/decoy_transcripts/gentrome.fa", **config),
            decoys=expand("{genome_dir}/{{assembly}}/decoy_transcripts/decoys.txt", **config),
        params:
            script=f"{config['rule_dir']}/../scripts/generateDecoyTranscriptome.sh",
        log:
            expand("{log_dir}/get_genome/{{assembly}}.decoy_transcripts.log", **config),
        message: explain_rule("decoy_transcripts")
        benchmark:
            expand("{benchmark_dir}/get_genome/{{assembly}}.decoy_transcripts.benchmark.txt", **config)[0]
        threads: 40
        resources:
            mem_gb=64,
        conda:
            "../envs/decoy.yaml"
        priority: 1
        shell:
            ("cpulimit --include-children -l {threads}00 -- " if config. get("cpulimit", True) else" ")+
            "sh {params.script} -j {threads} -g {input.genome} -a {input.gtf} -t {input.transcripts} -o $(dirname {output[0]}) > {log} 2>&1"

    rule salmon_decoy_aware_index:
        """
        Generate a decoy aware transcriptome index for Salmon.
        """
        input:
            gentrome=expand("{genome_dir}/{{assembly}}/decoy_transcripts/gentrome.fa", **config),
            decoys=expand("{genome_dir}/{{assembly}}/decoy_transcripts/decoys.txt", **config),
        output:
            directory(expand("{genome_dir}/{{assembly}}/index/{quantifier}_decoy_aware", **config)),
        log:
            expand("{log_dir}/{quantifier}_index/{{assembly}}.log", **config),
        benchmark:
            expand("{benchmark_dir}/{quantifier}_index/{{assembly}}.benchmark.txt", **config)[0]
        params:
            config["quantifier_index"],
        priority: 1
        threads: 10
        conda:
            "../envs/salmon.yaml"
        shell:
            """
            mkdir -p {output}
            
            salmon index -t {input.gentrome} -d {input.decoys} -i {output} {params} \
            --threads {threads} > {log} 2>&1
            """

    rule salmon_index:
        """
        Generate a transcriptome index for Salmon.
        """
        input:
            expand("{genome_dir}/{{assembly}}/{{assembly}}.transcripts.fa", **config),
        output:
            directory(expand("{genome_dir}/{{assembly}}/index/{quantifier}", **config)),
        log:
            expand("{log_dir}/{quantifier}_index/{{assembly}}.log", **config),
        benchmark:
            expand("{benchmark_dir}/{quantifier}_index/{{assembly}}.benchmark.txt", **config)[0]
        params:
            config["quantifier_index"],
        priority: 1
        threads: 10
        conda:
            "../envs/salmon.yaml"
        shell:
            """
            mkdir -p {output}
            
            salmon index -t {input} -i {output} {params} --threads {threads} > {log} 2>&1
            """

    rule salmon_quant:
        """
        Align reads against a transcriptome (index) with Salmon (mapping-based mode) and output a quantification file per sample.
        """
        input:
            reads=get_reads,
            index=get_salmon_index,
        output:
            dir=directory(expand("{result_dir}/{quantifier}/{{assembly}}-{{sample}}", **config)),
        log:
            expand("{log_dir}/{quantifier}_quant/{{assembly}}-{{sample}}.log", **config),
        message: explain_rule("salmon_quant")
        benchmark:
            expand("{benchmark_dir}/{quantifier}_quant/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
        params:
            input=(
                lambda wildcards, input: ["-r", input.reads]
                if sampledict[wildcards.sample]["layout"] == "SINGLE"
                else ["-1", input.reads[0], "-2", input.reads[1]]
            ),
            params=config["quantifier_flags"],
            reps=lambda wildcards, input: input  # help resolve changes in input files
        threads: 12
        resources:
            mem_gb=8,
        conda:
            "../envs/salmon.yaml"
        shell:
            """
            salmon quant -i {input.index} -l A {params.input} {params.params} -o {output.dir} \
            --threads {threads} > {log} 2>&1
            """


elif config["quantifier"] == "kallistobus":

    rule kallistobus_ref:
        """
        Make a genome index for kallistobus. This index is required for counting.
        """
        input:
            fa=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config),
            gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config),
        output:
            directory(expand("{genome_dir}/{{assembly}}/index/kallistobus/", **config)),
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
            basename=lambda wildcards, output: f"{output[0]}{wildcards.assembly}",
            options=config.get("ref")
        shell:
            """
            kb ref \
            {input.fa} {input.gtf} \
            -i {params.basename}.idx -g {params.basename}_t2g.txt -f1 {params.basename}_cdna.fa \
            -f2 {params.basename}_intron.fa \
            -c1 {params.basename}_cdna_t2c.txt -c2 {params.basename}_intron_t2c.txt \
            {params.options} > {log} 2>&1
            """


    def get_fastq_pair_reads(wildcards):
        """
        Extracts the correct combination of R1/R2 (trimmed and barcodes) for fastq_pair 
        based on Kallisto bustools settings.
        """
        reads = dict()
        assert sampledict[sample]["layout"] == "PAIRED", "Seq2science does not support scRNA-seq samples that are single-ended"
        read_id = get_bustools_rid(config.get("count"))
        #Determine mate for trimming
        if read_id == 0:
            reads["r1"] = expand("{trimmed_dir}/{{sample}}_{fqext1}_trimmed.{fqsuffix}.gz", **config)
            reads["r2"] = expand("{fastq_dir}/{{sample}}_{fqext2}.{fqsuffix}.gz", **config)
        elif read_id == 1:
            reads["r1"] = expand("{fastq_dir}/{{sample}}_{fqext1}.{fqsuffix}.gz", **config)
            reads["r2"] = expand("{trimmed_dir}/{{sample}}_{fqext2}_trimmed.{fqsuffix}.gz", **config)
        else:
            raise NotImplementedError
        return reads


    rule fastq_pair:    
        """
        fastq_pair re-writes paired-end fastq files to ensure that each read has a mate and 
        dsicards singleton reads. This step is required after scRNA trimming since we only trim the fastq 
        containing reads and not the barcode fastq. 
        """
        input:
            unpack(get_fastq_pair_reads)
        output:
            reads=temp(expand("{fastq_clean_dir}/{{sample}}_clean_{fqext}.{fqsuffix}.paired.fq", **config)),
            intermediates1=temp(expand("{fastq_clean_dir}/{{sample}}_clean_{fqext}.{fqsuffix}", **config)),
            intermediates2=temp(expand("{fastq_clean_dir}/{{sample}}_clean_{fqext}.{fqsuffix}.single.fq", **config))
        priority: 1
        log:
            expand("{log_dir}/fastq_pair/{{sample}}.log", **config),
        benchmark:
            expand("{benchmark_dir}/fastq_pair/{{sample}}.benchmark.txt", **config)[0]
        conda:
            "../envs/fastq-pair.yaml"
        params:
            options=config.get("fastq-pair",""),
            tused=lambda wildcards, input: "true" if "-t" in config.get("fastq-pair", "") else "false"
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
            """


    rule kallistobus_count:
        """
        Align reads against a transcriptome (index) with kallistobus and output a quantification file per sample.
        """
        input:
             barcodefile=config["barcodefile"],
             basedir=rules.kallistobus_ref.output,
             reads=rules.fastq_pair.output.reads
        output:
            dir=directory(expand("{result_dir}/{quantifier}/{{assembly}}-{{sample}}",**config))
        log:
            expand("{log_dir}/kallistobus_count/{{assembly}}-{{sample}}.log", **config),
        benchmark:
            expand("{benchmark_dir}/kallistobus_count/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
        priority: 1
        conda:
            "../envs/kallistobus.yaml"
        threads: 8
        message: explain_rule("kallistobus-count")
        resources:
            mem_gb=66,
        params:
            basename=lambda wildcards, input: f"{input.basedir[0]}/{wildcards.assembly}",
            options=config.get("count")
        shell:
            """
            kb count \
            -i {params.basename}.idx -w {input.barcodefile} \
            -t {threads} -g {params.basename}_t2g.txt \
            -o {output} -c1 {params.basename}_cdna_t2c.txt -c2 {params.basename}_intron_t2c.txt \
            {params.options} {input.reads} > {log} 2>&1
            """

    rule kb_seurat_pp:
        input:
            expand([f"{{result_dir}}/{{quantifier}}/{custom_assembly(treps.loc[trep, 'assembly'])}-{trep}" for trep in treps.index], **config)
        output:
            html=f"{config['result_dir']}/kb_seurat_pp/{{quantifier}}/{{assembly}}/kb_seurat_pp.html",
            qc_dir=directory(f"{config['result_dir']}/kb_seurat_pp/{{quantifier}}/{{assembly}}/qc")
        priority: 1
        conda:
            "../envs/kb_seurat_pp.yaml"
        params:
            isvelo=lambda wildcards, input: True if "--workflow" in config.get("count", "") else False
        message: explain_rule("kb_seurat_pp")    
        resources:
            R_scripts=1, # conda's R can have issues when starting multiple times
        script:
            f"{config['rule_dir']}/../scripts/knit_kb_seurat_pp.R"


elif config["quantifier"] == "htseq":

    rule htseq_count:
        """
        summarize reads to gene level. Outputs a counts table per bam file.
        """
        input:
            bam=expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config),
            gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config),
            required=_strandedness_report,
        output:
            expand("{counts_dir}/{{assembly}}-{{sample}}.counts.tsv", **config),
        params:
            strandedness=lambda wildcards: strandedness_to_quant(wildcards, "htseq"),
            user_flags=config["htseq_flags"]
        log:
            expand("{log_dir}/counts_matrix/{{assembly}}-{{sample}}.counts.log", **config),
        message: explain_rule("htseq_count")
        threads: 1
        conda:
            "../envs/gene_counts.yaml"
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
            bam=expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config),
            gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config),
            required=_strandedness_report,
        output:
            expand("{counts_dir}/{{assembly}}-{{sample}}.counts.tsv", **config),
        params:
            strandedness=lambda wildcards: strandedness_to_quant(wildcards, "featurecounts"),
            endedness=lambda wildcards: "" if sampledict[wildcards.sample]["layout"] == 'SINGLE' else "-p",
            user_flags=config["featurecounts_flags"],
        log:
            expand("{log_dir}/counts_matrix/{{assembly}}-{{sample}}.counts.log", **config),
        message: explain_rule("featurecounts_rna")
        threads: 1
        conda:
            "../envs/gene_counts.yaml"
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
             expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config),
        output:
             expand("{genome_dir}/{{assembly}}/{{assembly}}.DEXseq_annotation.gff", **config),
        log:
             expand("{log_dir}/counts_matrix/{{assembly}}.prepare_DEXseq_annotation.log", **config),
        conda:
             "../envs/dexseq.yaml"
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
            bam=expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config),
            bai=expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam.bai", **config),
            gff=expand("{genome_dir}/{{assembly}}/{{assembly}}.DEXseq_annotation.gff", **config),
            required=_strandedness_report,
        output:
            expand("{counts_dir}/{{assembly}}-{{sample}}.DEXSeq_counts.tsv", **config),
        log:
            expand("{log_dir}/counts_matrix/{{assembly}}-{{sample}}.DEXseq_counts.log", **config),
        message: explain_rule("dexseq")
        params:
            strandedness=lambda wildcards: strandedness_to_quant(wildcards, "dexseq"),
            endedness=lambda wildcards: "" if sampledict[wildcards.sample]["layout"] == 'SINGLE' else "-p yes",
        conda:
             "../envs/dexseq.yaml"
        shell:
             """
             current_conda_env=$(conda env list | grep \* | cut -d "*" -f2-)
             DEXseq_path=${{current_conda_env}}/lib/R/library/DEXSeq/python_scripts
             
             python ${{DEXseq_path}}/dexseq_count.py -f bam -r pos {params.endedness} -s {params.strandedness} {input.gff} {input.bam} {output} > {log} 2>&1
             """
