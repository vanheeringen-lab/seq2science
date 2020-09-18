import os
import shutil

onstart:
    # save a copy of the latest samples and config file(s) in the log_dir
    # skip this step on Jenkins, as it runs in parallel
    if os.getcwd() != config['log_dir'] and not os.getcwd().startswith('/var/lib/jenkins'):
        os.makedirs(config['log_dir'], exist_ok=True)
        for n, file in enumerate([config['samples']] + workflow.overwrite_configfiles):
            src = os.path.join(os.getcwd(), file)
            dst = os.path.join(config['log_dir'], os.path.basename(file) if n<2 else "profile.yaml")
            shutil.copy(src, dst)
onsuccess:
    if config.get("email") not in ["none@provided.com", "yourmail@here.com", None]:
        os.system(f"""echo "Succesful pipeline run! :)" | mail -s "The seq2science pipeline finished succesfully." {config["email"]} 2> /dev/null""")
onerror:
    if config.get("email") not in ["none@provided.com", "yourmail@here.com", None]:
        os.system(f"""echo "Unsuccessful pipeline run! :(" | mail -s "The seq2science pipeline finished prematurely..." {config["email"]} 2> /dev/null """)


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
keys_to_remove = ["fqext1", "fqext2", "macs2_types", "cpulimit",
                  "genome_types", "genomepy_temp", "bam_sort_mem",
                  ("biological_replicates", "condition" not in samples),
                  ("filter_bam_by_strand", "strandedness" not in samples),
                  ("technical_replicates", "replicates" not in samples),
                  ("tximeta", config.get("quantifier") != "salmon"),
                  ("deseq2", not config.get("contrasts")),
                  ("dge_dir", not config.get("contrasts")),
                  ("bigwig_dir", not config.get("create_trackhub")),
                  ("qc_dir", not config.get("create_qc_report"))]
keys = rmkeys(["samples"] + keys_to_remove, keys)
keys = ["samples"] + keys

for key in keys:
    if config[key] not in ["", False, 0, "None", "none@provided.com", "yourmail@here.com", "_custom"]:
        logger.info(f"{key: <23}: {config[key]}")

layouts = {sample: values['layout'] for sample, values in sampledict.items()}
logger.info(f"layout:                : {layouts}")


logger.info("\n\n")
