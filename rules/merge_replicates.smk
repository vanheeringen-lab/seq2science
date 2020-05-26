# dataframe with all technical replicates collapsed
cols = ['sample', 'assembly']
if 'replicate' in samples and config.get('technical_replicates')  != 'keep':
    cols = ['replicate', 'assembly']
if 'condition' in samples and config.get('biological_replicates') != 'keep':
    cols.append('condition')
if "descriptive_name" in samples:
    cols.append('descriptive_name')
if "control" in samples:
    cols.append("control")
treps = samples.reset_index()[cols].drop_duplicates().set_index(cols[0])

# dataframe with all replicates collapsed
breps = treps
if 'condition' in treps:
    breps = treps.reset_index(drop=True).drop_duplicates().set_index('condition')


# make a dict that returns the treps that belong to a brep
treps_from_brep = dict()
if "condition" in treps:
    for brep, row in breps.iterrows():
        assembly = row["assembly"]
        treps_from_brep[(brep, assembly)] = list(treps[(treps["assembly"] == assembly) & (treps["condition"] == brep)].index)
else:
    for brep, row in breps.iterrows():
        assembly = row["assembly"]
        treps_from_brep[(brep, assembly)] = [brep]

# and vice versa
brep_from_trep = dict()
for (brep, _assembly), _treps in treps_from_brep.items():
    brep_from_trep.update({trep: brep for trep in _treps})


def rep_to_descriptive(rep):
    """
    Return the descriptive name for a replicate.
    """
    if "descriptive_name" not in samples:
        return rep

    if rep in samples.index:
        return samples.loc[rep, "descriptive_name"]

    return rep


if 'replicate' in samples and config.get('technical_replicates') == 'merge':
    def get_merge_replicates(wildcards):
        output = dict()
        # since in the scATAC workflow an aditional rule called: "add_cell_id_2_fastq" adds the cell ID to the fastq and renames them 
        # to include the "_cell_ID_added" to the fastq file name, if the workflow that is ran is the scATAC one the input file names
        # of the merge replicate rule is changed to inlude the "_cell_ID_added" flagg. 
        scATAC = "_cell_ID_added" if get_workflow() == "scATAC" else ""
        output["reps"] = expand([f"{{trimmed_dir}}/{sample}{wildcards.fqext}_trimmed{scATAC}.{{fqsuffix}}.gz"
                         for sample in samples[samples['replicate'] == wildcards.replicate].index], **config)
        # make sure we make the fastqc report before moving our file
        if get_workflow() != "scATAC_seq" and len(output["reps"]) == 1 and config["create_qc_report"]:
            print('if get_workflow() != "scATAC_seq" doesnt work')
            output["qc"] = expand([f"{{qc_dir}}/fastqc/{sample}{wildcards.fqext}_trimmed{scATAC}_fastqc.zip"
                           for sample in samples[samples['replicate'] == wildcards.replicate].index], **config)
        return output

    rule merge_replicates:
        """
        Merge replicates (fastqs) simply by concatenating the files.
        
        Must happen after trimming due to trim-galore's automatic adapter trimming method 
        
        If a replicate has only 1 sample in it, rename and move instead.
        """
        input:
            unpack(get_merge_replicates)
        output:
            temp(sorted(expand("{trimmed_dir}/merged/{{replicate}}{{fqext}}_trimmed.{fqsuffix}.gz", **config)))
        wildcard_constraints:
            replicate=any_given('replicate'),
            fqext=f"_{config['fqext1']}|_{config['fqext2']}|" # nothing (SE), or fqext with an underscore (PE)
        log:
            expand("{log_dir}/merge_replicates/{{replicate}}{{fqext}}.log", **config)
        benchmark:
            expand("{benchmark_dir}/merge_replicates/{{replicate}}{{fqext}}.benchmark.txt", **config)[0]
        run:
            if len(input.reps) == 1:
                shell(f"mv {input.reps} {output} 2> {log}")
            else:
                for rep in input.reps:
                    rep_name = re.findall('\/([^\/_]+)_', rep)[-1]
                    print(rep, rep_name)
                    shell(f"""zcat {rep} | awk '{{{{gsub(/@/, "@{rep_name}:" )}}}}1' | gzip >> {output}""")
