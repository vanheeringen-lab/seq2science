if 'condition' in samples and config.get('combine_replicates', '') == 'merge':
    def get_merge_replicates(wildcards):
        return expand([f"{{trimmed_dir}}/{replicate}{wildcards.fqext}_trimmed.{{fqsuffix}}.gz"
               for replicate in samples[samples['condition'] == wildcards.condition].index], **config)

    rule merge_replicates:
        """
        Merge replicates (fastqs) simply by concatenating the files.
        """
        input:
            get_merge_replicates
        output:
            sorted(expand("{trimmed_dir}/merged/{{condition}}{{fqext}}_trimmed.{fqsuffix}.gz", **config))
        wildcard_constraints:
            fqext=".*",
            condition="[^/_]*"
        log:
            expand("{log_dir}/merge_replicates/{{condition}}{{fqext}}.log", **config)
        benchmark:
            expand("{benchmark_dir}/merge_replicates/{{condition}}{{fqext}}.benchmark.txt", **config)[0]
        shell:
            "cat {input} > {output} 2> {log}"
