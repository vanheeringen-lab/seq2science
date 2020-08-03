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
    assert breps.index.is_unique, "duplicate value found in breps"


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
        input_files = dict()
        input_files["reps"] = expand(
            [
                f"{{trimmed_dir}}/{sample}{wildcards.fqext}_trimmed.{{fqsuffix}}.gz"
                for sample in samples[samples["replicate"] == wildcards.replicate].index
            ],
            **config,
        )
        # make sure we make the fastqc report before moving our file
        if get_workflow() != "scATAC_seq" and len(input_files["reps"]) == 1 and config["create_qc_report"]:
            input_files["qc"] = expand(
                [
                    f"{{qc_dir}}/fastqc/{sample}{wildcards.fqext}_trimmed_fastqc.zip"
                    for sample in samples[samples["replicate"] == wildcards.replicate].index
                ],
                **config,
            )
        return input_files

    rule merge_replicates:
        """
        Merge replicates (fastqs) simply by concatenating the files. We also change the name of the read headers to 
        contain the name of the original replicate.

        Must happen after trimming due to trim-galore's automatic adapter trimming method 

        If a replicate has only 1 sample in it, simply move the file.
        """
        input:
            unpack(get_merge_replicates),
        output:
            temp(sorted(expand("{trimmed_dir}/merged/{{replicate}}{{fqext}}_trimmed.{fqsuffix}.gz", **config))),
        wildcard_constraints:
            replicate=any_given("replicate", "control"),
            fqext=f"_{config['fqext1']}|_{config['fqext2']}|", # nothing (SE), or fqext with an underscore (PE)
        log:
            expand("{log_dir}/merge_replicates/{{replicate}}{{fqext}}.log", **config),
        benchmark:
            expand("{benchmark_dir}/merge_replicates/{{replicate}}{{fqext}}.benchmark.txt", **config)[0]
        run:
            if len(input.reps) == 1:
                shell("mv {input.reps} {output} 2> {log}")
            else:
                for rep in input.reps:
                    rep_name = re.findall("\/([^\/_]+)_", rep)[-1]
                    shell(
                        """zcat {rep} | awk '{{if (NR%4==1) {{gsub(/^@/, "@{rep_name}:"); print}} else {{print}}}}' | gzip >> {output}"""
                    )
