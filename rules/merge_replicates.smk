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
        return expand([f"{{trimmed_dir}}/{sample}{wildcards.fqext}_trimmed.{{fqsuffix}}.gz"
               for sample in samples[samples['replicate'] == wildcards.replicate].index], **config)

    rule merge_replicates:
        """
        Merge replicates (fastqs) simply by concatenating the files.
        
        Must happen after trimming due to trim-galore's automatic adapter trimming method 
        
        If a replicate has only 1 sample in it, symlink it instead.
        """
        input:
            get_merge_replicates
        output:
            temp(sorted(expand("{trimmed_dir}/merged/{{replicate}}{{fqext}}_trimmed.{fqsuffix}.gz", **config)))
        wildcard_constraints:
            replicate=any_given('replicate'),
            fqext=f"_{config['fqext1']}|_{config['fqext2']}|" # nothing (SE), or fqext with an underscore (PE)
        log:
            expand("{log_dir}/merge_replicates/{{replicate}}{{fqext}}.log", **config)
        benchmark:
            expand("{benchmark_dir}/merge_replicates/{{replicate}}{{fqext}}.benchmark.txt", **config)[0]
        shell:
            """
            arr=({input})
            if [ ${{#arr[@]}} -eq 1 ]; then
                echo '\nlinking file:\n{input}' > {log}
                ln -s {input} {output}  2> {log}
            else 
                echo '\nconcatenating files:\n{input}' > {log}
                cat {input} > {output} 2> {log}
            fi
            """
