# dataframe with all technical replicates collapsed
cols = ["sample", "assembly"]
subset = ["sample", "assembly"]
if "replicate" in samples:
    cols = ["replicate", "assembly"]
    subset = ["replicate", "assembly"]
if "condition" in samples:
    cols.append("condition")
    subset.append("condition")
if "control" in samples:
    cols.append("control")
if "colors" in samples:
    cols.append("colors")

treps = samples.reset_index()[cols].drop_duplicates(subset=subset).set_index(cols[0])
assert treps.index.is_unique, "duplicate value found in treps"

# treps that came from a merge
merged_treps = [trep for trep in treps.index if trep not in samples.index]
merged_treps_single = [trep for trep in merged_treps if sampledict[trep]["layout"] == "SINGLE"]
merged_treps_paired = [trep for trep in merged_treps if sampledict[trep]["layout"] == "PAIRED"]

# all samples (including controls)
all_samples = [sample for sample in samples.index]
if "control" in samples.columns:
    all_samples.extend(samples["control"].dropna().tolist())

# dataframe with all replicates collapsed
breps = treps
if "condition" in treps and config.get("biological_replicates", "keep") != "keep":
    breps = treps.reset_index(drop=True).drop_duplicates(subset=subset[1:]).set_index("condition")


# make a dict that returns the treps that belong to a brep
treps_from_brep = dict()
if "condition" in treps and config.get("biological_replicates", "keep") != "keep":
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
            rep = samples[samples.condition == rep].condition[0]
        else:
            if "replicate" in samples:
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
