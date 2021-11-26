"""
all rules/logic related to the merging of technical replicates (and a bit for biological replicates) should be here.
"""

# dataframe with all technical replicates collapsed
cols = ["sample", "assembly"]
subset = ["sample", "assembly"]
if "technical_replicates" in samples:
    cols = ["technical_replicates", "assembly"]
    subset = ["technical_replicates", "assembly"]
if "biological_replicates" in samples:
    cols.append("biological_replicates")
    subset.append("biological_replicates")
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
if "biological_replicates" in treps:
    breps = treps.reset_index(drop=True).drop_duplicates(subset=subset[1:]).set_index("biological_replicates")


# make a dict that returns the treps that belong to a brep
treps_from_brep = dict()
if "biological_replicates" in treps:
    for brep, row in breps.iterrows():
        assembly = row["assembly"]
        treps_from_brep[(brep, assembly)] = list(
            treps[(treps["assembly"] == assembly) & (treps["biological_replicates"] == brep)].index
        )
        treps_from_brep[(brep, assembly + config.get("custom_assembly_suffix", ""))] = list(
            treps[(treps["assembly"] == assembly) & (treps["biological_replicates"] == brep)].index
        )
else:
    for brep, row in breps.iterrows():
        assembly = row["assembly"]
        treps_from_brep[(brep, assembly)] = [brep]
        treps_from_brep[(brep, assembly + config.get("custom_assembly_suffix", ""))] = [brep]

# and vice versa
brep_from_trep = dict()
for (brep, _assembly), _treps in treps_from_brep.items():
    brep_from_trep.update({trep: brep for trep in _treps})


def rep_to_descriptive(rep, brep=False):
    """
    Return the descriptive name for a replicate.
    """
    if "descriptive_name" in samples:
        if brep and "biological_replicates" in samples:
            rep = samples[samples.biological_replicates == rep].biological_replicates[0]
        else:
            if "technical_replicates" in samples:
                col = samples.technical_replicates
            else:
                col = samples.index
            rep = samples[col == rep].descriptive_name[0]
    return rep


if "technical_replicates" in samples:

    def get_merge_replicates(wildcards):
        input_files = list()
        for sample in samples[samples["technical_replicates"] == wildcards.replicate].index:
            if get_workflow() == "scrna_seq":
                input_files.append(f"{{fastq_clean_dir}}/{sample}_clean{wildcards.fqext}.{{fqsuffix}}.paired.fq")
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
            fqext=f"_{config['fqext1']}|_{config['fqext2']}|", # nothing (SE), or fqext with an underscore (PE)
        params:
            reps=lambda wildcards, input: input  # help resolve changes in input files
        log:
            expand("{log_dir}/merge_replicates/{{replicate}}{{fqext}}.log", **config),
        benchmark:
            expand("{benchmark_dir}/merge_replicates/{{replicate}}{{fqext}}.benchmark.txt", **config)[0]
        run:
            # all workflows use gzipped fastq files, however scRNA-seq uses unzipped fastqs
            # due to trimming + pairing, that's why we set the printcmd
            printcmd = "cat" if get_workflow() == "scrna_seq" else "zcat"
            for rep in input:
                rep_name = re.findall("\/([^\/_]+)_", rep)[-1]

                shell(
                    """{printcmd} {rep} | awk '{{if (NR%4==1) {{gsub(/^@/, "@{rep_name}:"); print}} else {{print}}}}' | gzip >> {output}"""
                )
