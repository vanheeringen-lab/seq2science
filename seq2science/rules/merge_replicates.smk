"""
all rules/logic related to the merging of technical replicates (and a bit for biological replicates) should be here.
"""
if "technical_replicates" in samples:
    import re


    def get_merge_replicates(wildcards):
        input_files = list()
        for sample in samples[samples["technical_replicates"] == wildcards.replicate].index:
            if WORKFLOW == "scrna_seq":
                input_files.append(f"{{fastq_clean_dir}}/{sample}_clean{wildcards.fqext}.{{fqsuffix}}.paired.fq.gz")
            else:
                input_files.append(f"{{trimmed_dir}}/{sample}{wildcards.fqext}_trimmed.{{fqsuffix}}.gz")

        input_files = expand(input_files, **config)
        if len(input_files) == 0:
            return ["Something went wrong, and we tried to merge replicates but there were no replicates?!"]

        return input_files

    # scrna-seq uses the fastq paired files
    if config["trimmer"] == "fastp":

        ruleorder: merge_replicates > fastp_PE > fastp_SE


    elif config["trimmer"] == "trimgalore":

        ruleorder: merge_replicates > trimgalore_PE > trimgalore_SE

    # true treps are treps combined of 2 samples or more
    true_treps = [trep for trep in treps.index if trep not in samples.index]

    rule merge_replicates:
        """
        Merge replicates (fastqs) simply by concatenating the files. We also change the name of the read headers to 
        contain the name of the original replicate.

        Must happen after trimming due to trim-galore's automatic adapter trimming method 

        If a replicate has only 1 sample in it, simply move the file.
        """
        input:
            get_merge_replicates,
        output:
            temp(sorted(expand("{trimmed_dir}/{{replicate}}{{fqext}}_trimmed.{fqsuffix}.gz", **config))),
        wildcard_constraints:
            replicate="|".join(true_treps) if len(true_treps) else "$a",
            fqext=f"_{config['fqext1']}|_{config['fqext2']}|",  # nothing (SE), or fqext with an underscore (PE)
        params:
            reps=lambda wildcards, input: input,  # help resolve changes in input files
        log:
            expand("{log_dir}/merge_replicates/{{replicate}}{{fqext}}.log", **config),
        benchmark:
            expand("{benchmark_dir}/merge_replicates/{{replicate}}{{fqext}}.benchmark.txt", **config)[0]
        run:
            for rep in input:
                rep_name = re.findall("/([^/_]+)_", rep)[-1]

                shell(
                    """zcat {rep} | awk '{{if (NR%4==1) {{gsub(/^@/, "@{rep_name}:"); print}} else {{print}}}}' | gzip >> {output}"""
                )
