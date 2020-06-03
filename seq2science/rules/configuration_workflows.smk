# apply workflow specific changes...
# ...for atac-seq
if config.get('peak_caller', False):
    config['peak_caller'] = {k: v for k,v in config['peak_caller'].items()}

    # if genrich is peak caller, make sure to not double shift reads
    if 'genrich' in config['peak_caller']:
        # always turn of genrich shift, since we handle that with deeptools
        if '-j' in config['peak_caller']['genrich'] and not '-D' in config['peak_caller']['genrich']:
            config['peak_caller']['genrich'] += ' -D'

    # if hmmratac peak caller, check if all samples are paired-end
    if 'hmmratac' in config['peak_caller']:
        assert all([config['layout'][sample] == 'PAIRED' for sample in samples.index]), \
        "HMMRATAC requires all samples to be paired end"

    config['macs2_types'] = ['control_lambda.bdg', 'peaks.xls', 'treat_pileup.bdg']
    if 'macs2' in config['peak_caller']:
        params = config['peak_caller']["macs2"].split(" ")
        invalid_params = ["-t", "--treatment", "-c", "--control", "-n", "--name", "--outdir", "-f",
                          "--format", "-g", "--gsize", "-p", "--pvalue"]
        assert not any(val in params for val in invalid_params), f"You filled in a parameter for macs2 which the " \
                                                                 f"pipeline does not support. Unsupported params are:" \
                                                                 f"{invalid_params}."

        config["macs_cmbreps"] = ""
        cmbreps_params = ["-q", "--qvalue", "--min-length", "--max-gap", "--broad-cutoff"]
        for param in cmbreps_params:
            if param in params:
                idx = params.index(param) + 1
                if param == "-q" or param == "--qvalue":
                    val = -math.log(float(params[idx]), 10)
                    config["macs_cmbreps"] += f" -c {val} "
                else:
                    config["macs_cmbreps"] += f" {param} {params[idx]} "

        if '--broad' in config['peak_caller']['macs2']:
            config['macs2_types'].extend(['peaks.broadPeak', 'peaks.gappedPeak'])
        else:
            config['macs2_types'].extend(['summits.bed', 'peaks.narrowPeak'])


# ...for alignment and rna-seq
for conf_dict in ['aligner', 'quantifier', 'diffexp']:
    if config.get(conf_dict, False):
        dict_key = list(config[conf_dict].keys())[0]
        for k, v in list(config[conf_dict].values())[0].items():
            config[k] = v
        config[conf_dict] = dict_key

# ...for alignment
if config.get('bam_sorter', False):
    config['bam_sort_order'] = list(config['bam_sorter'].values())[0]
    config['bam_sorter'] = list(config['bam_sorter'].keys())[0]


# make sure that our samples.tsv and configuration work together...
# ...on biological replicates
if 'condition' in samples and config.get('biological_replicates', '') != 'keep':
    if 'hmmratac' in config['peak_caller']:
        assert config.get('biological_replicates', '') == 'idr', \
        f'HMMRATAC peaks can only be combined through idr'

    for condition in set(samples['condition']):
        for assembly in set(samples[samples['condition'] == condition]['assembly']):
            if 'replicate' in samples and config.get('technical_replicates') == 'merge':
                nr_samples = len(set(samples[(samples['condition'] == condition) & (samples['assembly'] == assembly)]['replicate']))
            else:
                nr_samples = len(samples[(samples['condition'] == condition) & (samples['assembly'] == assembly)])

            if config.get('biological_replicates', '') == 'idr':
                assert nr_samples <= 2,\
                f'For IDR to work you need two samples per condition, however you gave {nr_samples} samples for'\
                f' condition {condition} and assembly {assembly}'

# ...on DE contrasts
def parse_DE_contrasts(de_contrast):
    """
    Extract batch and contrast groups from a DE contrast design
    """
    original_contrast = de_contrast

    # remove whitespaces (and '~'s if used)
    de_contrast = de_contrast.replace(" ", "").replace("~", "")

    # split and store batch effect
    batch = None
    if '+' in de_contrast:
        batch =  de_contrast.split('+')[0]
        de_contrast = de_contrast.split('+')[1]

    # parse contrast
    parsed_contrast = de_contrast.split('_')
    return original_contrast, parsed_contrast, batch

if config.get('contrasts', False):
    # check differential gene expression contrasts
    old_contrasts = list(config["contrasts"])
    for contrast in old_contrasts:
        original_contrast, parsed_contrast, batch = parse_DE_contrasts(contrast)

        # Check if the column names can be recognized in the contrast
        assert parsed_contrast[0] in samples.columns and parsed_contrast[0] not in ["sample", "assembly"], \
            (f'\nIn contrast design "{original_contrast}", "{parsed_contrast[0]} ' +
             f'does not match any valid column name in {config["samples"]}.\n')
        if batch is not None:
            assert batch in samples.columns and batch not in ["sample", "assembly"], \
                (f'\nIn contrast design "{original_contrast}", the batch effect "{batch}" ' +
                 f'does not match any valid column name in {config["samples"]}.\n')

        # Check if the groups described by the contrast can be identified and found in samples.tsv
        l = len(parsed_contrast)
        assert l < 4, ("\nA differential expression contrast couldn't be parsed correctly.\n" +
                       f"{str(l-1)} groups were found in '{original_contrast}' " +
                       f"(groups: {', '.join(parsed_contrast[1:])}).\n\n" +
                       f'Tip: do not use underscores in the columns of {config["samples"]} referenced by your contrast.\n')
        if l == 1:
            # check if contrast column has exactly 2 factor levels (per assembly)
            tmp = samples[['assembly', parsed_contrast[0]]].dropna()
            factors = pd.DataFrame(tmp.groupby('assembly')[parsed_contrast[0]].nunique())
            assert all(factors[parsed_contrast[0]] == 2),\
                ('\nYour contrast design, ' + original_contrast +
                 ', contains only a column name (' + parsed_contrast[0] +
                 '). \nIf you wish to compare all groups in this column, add a reference group. ' +
                 'Number of groups found (per assembly): \n' + str(factors[parsed_contrast[0]]))
        else:
            # check if contrast column contains the groups
            for group in parsed_contrast[1:]:
                if group != 'all':
                    assert str(group) in [str(i) for i in samples[parsed_contrast[0]].tolist()],\
                        ('\nYour contrast design contains group ' + group +
                        ' which cannot be found in column ' + parsed_contrast[0] +
                         ' of ' + config["samples"] + '.\n')


def sieve_bam(configdict):
    """
    helper function to check whether or not we use rule sieve_bam
    """
    return configdict.get('min_mapping_quality', 0) > 0 or \
           configdict.get('tn5_shift', False) or \
           configdict.get('remove_blacklist', False) or \
           configdict.get('remove_mito', False)
