# dataframe with all technical replicates collapsed
cols = ["sample", "assembly"]
if "replicate" in samples:
    cols = ["replicate", "assembly"]
if "condition" in samples:
    cols.append("condition")
if "control" in samples:
    cols.append("control")
treps = samples.reset_index()[cols].drop_duplicates().set_index(cols[0])
assert treps.index.is_unique, "duplicate value found in treps"

# dataframe with all replicates collapsed
breps = treps
if "condition" in treps:
    breps = treps.reset_index(drop=True).drop_duplicates().set_index("condition")


# make a dict that returns the treps that belong to a brep
treps_from_brep = dict()
if "condition" in treps:
    for brep, row in breps.iterrows():
        assembly = row["assembly"]
        treps_from_brep[(brep, assembly)] = list(
            treps[(treps["assembly"] == assembly) & (treps["condition"] == brep)].index
        )
else:
    for brep, row in breps.iterrows():
        assembly = row["assembly"]
        treps_from_brep[(brep, assembly)] = [brep]

# and vice versa
brep_from_trep = dict()
for (brep, _assembly), _treps in treps_from_brep.items():
    brep_from_trep.update({trep: brep for trep in _treps})


def rep_to_descriptive(rep, brep=False):
    """
    Return the descriptive name for a replicate.
    """
    if "descriptive_name" in samples:
        if brep and "condition" in samples:
            col = samples.condition
        elif "replicate" in samples:
            col = samples.replicate
        else:
            col = samples.index
        rep = samples[col == rep].descriptive_name[0]
    return rep


if "replicate" in samples:

    def get_merge_replicates(wildcards):
        input_files = expand(
            [
                f"{{trimmed_dir}}/{sample}{wildcards.fqext}_trimmed.{{fqsuffix}}.gz"
                for sample in samples[samples["replicate"] == wildcards.replicate].index
            ],
            **config,
        )
        if len(input_files) == 0:
            return ["Something went wrong, and we tried to merge replicates but there were no replicates?!"]

        return input_files

    if config["trimmer"] == "fastp":
        ruleorder: merge_replicates > fastp_PE > fastp_SE
    elif config["trimmer"] == "trimgalore":
        ruleorder: merge_replicates > trimgalore_PE > trimgalore_SE

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
            replicate=any_given("replicate"),
            fqext=f"_{config['fqext1']}|_{config['fqext2']}|", # nothing (SE), or fqext with an underscore (PE)
        log:
            expand("{log_dir}/merge_replicates/{{replicate}}{{fqext}}.log", **config),
        benchmark:
            expand("{benchmark_dir}/merge_replicates/{{replicate}}{{fqext}}.benchmark.txt", **config)[0]
        run:
            for rep in input:
                rep_name = re.findall("\/([^\/_]+)_", rep)[-1]
                shell(
                    """zcat {rep} | awk '{{if (NR%4==1) {{gsub(/^@/, "@{rep_name}:"); print}} else {{print}}}}' | gzip >> {output}"""
                )
