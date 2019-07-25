rule cutadapt:
    input:
        dir=download_fastq(expand("{result_dir}/FASTQ/{{sample}}/", **config)),
        adapt=expand("{result_dir}/FASTQ/{{sample}}/adapter.txt", **config)
    output:
        dir=directory(expand("{result_dir}/trimmed/{{sample}}", **config)),
        qc=expand("{result_dir}/trimmed_qc/{{sample}}/{{sample}}_trimming_report.txt", **config)
    log:
        expand("{log_dir}/cutadapt/{{sample}}.log", **config)
    params:
    # TODO: check in config if read1/read2 is sorted
        input=sorted(expand("{result_dir}/FASTQ/{{sample}}/*.{fqsuffix}.gz", **config)),
        read2=expand("{result_dir}/FASTQ/{{sample}}/{{sample}}_{fqext2}.{fqsuffix}.gz", **config)[0],
    # TODO: add config params
        config="-a AGATCGGAAGAGCAGATCGGAAGAGCAGATCGGAAGAGC"
    conda:
        "../envs/trim_auto.yaml"
    threads: 8
    shell:
        """
        read1=$(echo {params.input} | cut -f1 -d' ')
        output="-o $read1"

        # if paired end then add read2 to the output
        if [[ -f {params.read2} ]]; then
            output="$output -p {params.read2}"
        fi;

        # make sure that the output gets written to the trimmed folder, not the fastq folder
        output="${{output//FASTQ/trimmed}}"

        # remove directory on premature end
        outdir=$(dirname "${{read1//FASTQ/trimmed}}")
        trap \"rm -f $outdir\" INT;
        # now cut adapters
        mkdir $outdir &
        cpulimit --include-children -l {threads}00 --\
        cutadapt -j {threads} {params.config} $output {params.input} 1> {output.qc} 2> {log}
        """


# FIXME: this obviously does not work!
# TODO: pipe? Or temp output?
rule detect_adapter:
    input:
        download_fastq(expand("{result_dir}/FASTQ/{{sample}}/", **config))
    output:
        expand("{result_dir}/FASTQ/{{sample}}/adapter.txt", **config)
    log: # TODO
        expand("{log_dir}/detect_adapter/{{sample}}.log", **config)
    run:
        with open(output[0], "w") as text_file:
            text_file.write("AGATCGGAAGAGC")
