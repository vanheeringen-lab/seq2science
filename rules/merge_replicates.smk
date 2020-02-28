# def replicates_or_samples_per_assembly():
#     """returns a dict of samples/replicates (values) per assembly (key)"""
#     rspa = {}
#     for assembly in set(samples['assembly']):
#         rspa[assembly] = []
#
#         if 'replicate' in samples and config.get('technical_replicates') == 'merge':
#             for fname in set(samples[samples['assembly'] == assembly].replicate):
#                 rspa[assembly].append(fname)
#         else:
#             for fname in set(samples[samples['assembly'] == assembly].index):
#                 rspa[assembly].append(fname)
#     return rspa
#
#
# def peakfiles_per_assembly():
#     ppa = replicates_or_samples_per_assembly()
#     if 'condition' in samples and config.get('biological_replicates') != 'keep':
#         s2 = samples[]
#
#         for fname in set(samples[samples['assembly'] == assembly].condition):
#             rspa[assembly].append(fname)
#     else:
#         return ppa
#
# print(peakfiles_per_assembly())
# exit(0)

# dataframe with all technical replicates collapsed
cols = ['sample', 'assembly']
if 'replicate' in samples: #and config.get('technical_replicates') == 'merge':
    cols = ['replicate', 'assembly']
if 'condition' in samples: #and config.get('biological_replicates') != 'keep':
    cols.append('condition')
treps = samples.reset_index()[cols].drop_duplicates().set_index(cols[0])

# dataframe with all replicates collapsed
breps = treps
if 'condition' in samples: # and config.get('biological_replicates') != 'keep':
    breps = treps.reset_index(drop=True).drop_duplicates().set_index('condition')

# print(treps)
# print(breps)
# exit(0)

# # dataframe with all technical replicates collapsed
# treps = samples
# if 'replicate' in samples and config.get('technical_replicates') == 'merge':
#     treps = samples.reset_index(drop=True).drop_duplicates().set_index('replicate')
#
# # dataframe with all replicates collapsed
# breps = treps
# if 'condition' in samples and config.get('biological_replicates') != 'keep':
#     breps = treps.reset_index(drop=True).drop_duplicates().set_index('condition')
#
# print(treps)
# print(breps)

if 'replicate' in samples and config.get('technical_replicates') == 'merge':
    def get_merge_replicates(wildcards):
        return expand([f"{{trimmed_dir}}/{sample}{wildcards.fqext}_trimmed.{{fqsuffix}}.gz"
               for sample in samples[samples['replicate'] == wildcards.replicate].index], **config)

    rule merge_replicates:
        """
        Merge replicates (fastqs) simply by concatenating the files.
        
        If a replicate has only 1 sample in it, symlink it instead.
        """
        input:
            get_merge_replicates
        output:
            sorted(expand("{trimmed_dir}/merged/{{replicate}}{{fqext}}_trimmed.{fqsuffix}.gz", **config))
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
                ln -s {input} {output}  2> {log}
            else 
                cat {input} > {output} 2> {log}
            fi
            """
