from snakemake.logging import logger


# apply workflow specific changes...
# ...for atac-seq
if config.get("peak_caller", False):
    config["peak_caller"] = {k: v for k, v in config["peak_caller"].items()}

    # if genrich is peak caller, make sure to not double shift reads
    if "genrich" in config["peak_caller"]:
        # always turn of genrich shift, since we handle that with deeptools
        if "-j" in config["peak_caller"]["genrich"] and not "-D" in config["peak_caller"]["genrich"]:
            config["peak_caller"]["genrich"] += " -D"

    # if hmmratac peak caller, check if all samples are paired-end
    if "hmmratac" in config["peak_caller"]:
        assert all(
            [config["layout"][sample] == "PAIRED" for sample in samples.index]
        ), "HMMRATAC requires all samples to be paired end"

    config["macs2_types"] = ["control_lambda.bdg", "peaks.xls", "treat_pileup.bdg"]
    if "macs2" in config["peak_caller"]:
        params = config["peak_caller"]["macs2"].split(" ")
        invalid_params = [
            "-t",
            "--treatment",
            "-c",
            "--control",
            "-n",
            "--name",
            "--outdir",
            "-f",
            "--format",
            "-g",
            "--gsize",
            "-p",
            "--pvalue",
        ]
        assert not any(val in params for val in invalid_params), (
            f"You filled in a parameter for macs2 which the "
            f"pipeline does not support. Unsupported params are:"
            f"{invalid_params}."
        )

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

        if "--broad" in config["peak_caller"]["macs2"]:
            config["macs2_types"].extend(["peaks.broadPeak", "peaks.gappedPeak"])
        else:
            config["macs2_types"].extend(["summits.bed", "peaks.narrowPeak"])


# ...for alignment and rna-seq
for conf_dict in ["aligner", "quantifier", "diffexp"]:
    if config.get(conf_dict, False):
        dict_key = list(config[conf_dict].keys())[0]
        for k, v in list(config[conf_dict].values())[0].items():
            config[k] = v
        config[conf_dict] = dict_key

# ...for alignment
if config.get("bam_sorter", False):
    config["bam_sort_order"] = list(config["bam_sorter"].values())[0]
    config["bam_sorter"] = list(config["bam_sorter"].keys())[0]


# make sure that our samples.tsv and configuration work together...
# ...on biological replicates
if "condition" in samples:
    if "hmmratac" in config["peak_caller"]:
        assert config.get("biological_replicates", "") == "idr", f"HMMRATAC peaks can only be combined through idr"

    for condition in set(samples["condition"]):
        for assembly in set(samples[samples["condition"] == condition]["assembly"]):
            if "replicate" in samples:
                nr_samples = len(
                    set(samples[(samples["condition"] == condition) & (samples["assembly"] == assembly)]["replicate"])
                )
            else:
                nr_samples = len(samples[(samples["condition"] == condition) & (samples["assembly"] == assembly)])

            if config.get("biological_replicates", "") == "idr":
                assert nr_samples <= 2, (
                    f"For IDR to work you need two samples per condition, however you gave {nr_samples} samples for"
                    f" condition {condition} and assembly {assembly}"
                )

# ...on DE contrasts
if config.get("contrasts"):

    def parse_de_contrasts(de_contrast):
        """
        Extract contrast column, groups and batch effect column from a DE contrast design

        Accepts a string containing a DESeq2 contrast design

        Returns
        parsed_contrast, lst: the contrast components
        batch, str: the batch effect column or None
        """
        # remove whitespaces (and '~'s if used)
        de_contrast = de_contrast.replace(" ", "").replace("~", "")

        # split and store batch effect
        batch = None
        if "+" in de_contrast:
            batch, de_contrast = de_contrast.split("+")[0:2]

        # parse contrast column and groups
        parsed_contrast = de_contrast.split("_")
        return parsed_contrast, batch

    # check differential gene expression contrasts
    for contrast in list(config["contrasts"]):
        parsed_contrast, batch = parse_de_contrasts(contrast)
        column_name = parsed_contrast[0]

        # Check if the column names can be recognized in the contrast
        assert column_name in samples.columns and column_name not in ["sample", "assembly"], (
            f'\nIn contrast design "{contrast}", "{column_name} '
            + f'does not match any valid column name in {config["samples"]}.\n'
        )
        if batch is not None:
            assert batch in samples.columns and batch not in ["sample", "assembly"], (
                f'\nIn contrast design "{contrast}", the batch effect "{batch}" '
                + f'does not match any valid column name in {config["samples"]}.\n'
            )

        # Check if the contrast contains the number of components we can parse
        components = len(parsed_contrast)
        assert components > 1, (
            f"\nA differential expression contrast is too short: '{contrast}'.\n"
            "Specify at least one reference group to compare against.\n"
        )
        assert components < 4, (
            "\nA differential expression contrast couldn't be parsed correctly.\n"
            + f"{str(components-1)} groups were found in '{contrast}' "
            + f"(groups: {', '.join(parsed_contrast[1:])}).\n\n"
            + f'Tip: do not use underscores in the columns of {config["samples"]} referenced by your contrast.\n'
        )

        # Check if the groups in the contrast can be found in the right column of the samples.tsv
        for group in parsed_contrast[1:]:
            if group != "all":
                assert str(group) in [str(i) for i in samples[column_name].tolist()], (
                    f"\nYour contrast design contains group {group} which cannot be found "
                    f'in column {column_name} of {config["samples"]}.\n'
                )


def sieve_bam(configdict):
    """
    helper function to check whether or not we use rule sieve_bam
    """
    return (
        configdict.get("min_mapping_quality", 0) > 0
        or configdict.get("tn5_shift", False)
        or configdict.get("remove_blacklist", False)
        or configdict.get("remove_mito", False)
    )


def rmkeys(del_list, target_list):
    """
    remove all elements in del_list from target_list
    each element may be a tuple with an added condition
    """
    for element in del_list:
        if isinstance(element, str) and element in target_list:
            target_list.remove(element)
        elif element[1] and element[0] in target_list:
            target_list.remove(element[0])
    return target_list


# after all is done, log (print) the configuration
logger.info("CONFIGURATION VARIABLES:")

# sort config: samples.tsv & directories first, alphabetized second
keys = sorted(config.keys())
dir_keys = []
other_keys = []
for key in keys:
    if key.endswith("_dir"):
        dir_keys.append(key)
    else:
        other_keys.append(key)
keys = dir_keys + other_keys

# remove superfluous keys
keys_to_remove = ["samples", "layout", "fqext1", "fqext2", "macs2_types",
                  "cpulimit", "genome_types", "genomepy_temp", "bam_sort_mem",
                  ("biological_replicates", "condition" not in samples),
                  ("filter_bam_by_strand", "strandedness" not in samples),
                  ("technical_replicates", "replicates" not in samples),
                  ("tximeta", config.get("quantifier") is not "salmon")]
keys = rmkeys(keys_to_remove, keys)
keys = ["samples"] + keys + ["layout"]

for key in keys:
    if config[key] not in ["", False, 0, "None", "none@provided.com"]:
        logger.info(f"{key: <23}: {config[key]}")
logger.info("\n\n")
