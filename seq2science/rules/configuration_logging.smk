"""
all logic related to logging to stdout/file should be here.
"""

import os
import shutil
import textwrap

from tabulate import tabulate
from snakemake.logging import logger


onstart:
    # save a copy of the latest samples and config file(s) in the log_dir
    # skip this step on Jenkins, as it runs in parallel
    if os.getcwd() != config["log_dir"] and not os.getcwd().startswith("/var/lib/jenkins"):
        os.makedirs(config["log_dir"], exist_ok=True)
        for n, file in enumerate([config["samples"]] + workflow.overwrite_configfiles):
            src = os.path.join(os.getcwd(), file)
            dst = os.path.join(config["log_dir"], os.path.basename(file) if n < 2 else "profile.yaml")
            shutil.copy(src, dst)


onsuccess:
    # email that seq2science finished
    if config.get("email") not in ["none@provided.com", "yourmail@here.com", None]:
        os.system(
            f"""echo "Succesful pipeline run! :)" | mail -s "The seq2science pipeline finished succesfully." {config["email"]} 2> /dev/null"""
        )

    # be happy that it finished succesfully
    # and indicate where important files are located
    message = f"Nice, a succesful run! Check out the docs for help with the results: https://vanheeringen-lab.github.io/seq2science/content/workflows/{WORKFLOW}.html. "
    if WORKFLOW != "download_fastq":
        if config.get("create_qc_report"):
            rand_assembly = list(ORI_ASSEMBLIES.keys())[0]
            message += f"Make sure to check out the QC report, it can be found at {config['qc_dir']}/multiqc_{rand_assembly}.html. "
        if WORKFLOW not in ("scatac-seq", "scrna-seq") and config.get("create_trackhub") and "trackhub_dir" in config:
            message += f"And make sure to publicly host the trackhub for UCSC. The hub is located in {config['trackhub_dir']}/."

    message = textwrap.wrap(message, 70, break_long_words=False, break_on_hyphens=False)
    dancer = """
      ⊂_ヽ
      　 ＼＼
      　　 ＼( ͡° ͜ʖ ͡°)  --  
       　　　 >　⌒ヽ       
       　　　/ 　 へ＼      
       　　 /　　/　＼＼    
       　　 ﾚ　ノ　　 ヽ_つ  
       　　/　/            
       　 /　/|            
       　(　(ヽ            
       　|　|、＼
       　| 丿 ＼ ⌒)
       　| |　　) /
        ノ )　　Lﾉ
       (_／
    """

    # wrap the dancer and the message
    dance_message = ""
    offset = 3
    for i, msg in enumerate(dancer.split("\n")):
        dance_message += msg
        if i >= offset and (i - offset) < len(message):
            dance_message += message[i - offset]
        dance_message += "\n"
    logger.info(dance_message)

onerror:
    logger.info(
        """


      __    __    __  ____  ____  _   
     /  \  /  \  /  \(  _ \/ ___)/ \  
    (  O )(  O )(  O )) __/\___ \\\\_/  
     \__/  \__/  \__/(__)  (____/(_)  


    One or more rules did not finish as expected!


    Please take a look at the log files of the failed rule(s), and our Frequently Asked Questions: 
    https://vanheeringen-lab.github.io/seq2science/content/faq.html

    If that does not help you, don't be afraid to reach out to us. 
    The easiest way would be to make an issue on our github page: 
    https://github.com/vanheeringen-lab/seq2science/issues

    """
    )
    if config.get("email") not in ["none@provided.com", "yourmail@here.com", None]:
        os.system(
            f"""echo "Unsuccessful pipeline run! :(" | mail -s "The seq2science pipeline finished prematurely..." {config["email"]} 2> /dev/null """
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


if not config.get("no_config_log"):
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

    # check if aligners are used at all
    no_aligners = config.get("quantifier") in ["salmon", "kallistobus"] and not config.get("create_trackhub")

    # remove superfluous keys
    keys_to_remove = [
        "fqext1",
        "fqext2",
        "macs2_types",
        "cpulimit",
        "genome_types",
        "genomepy_temp",
        "bam_sort_mem",
        "benchmark_dir",
        "sra_dir",
        "trimmed_dir",
        "result_dir",
        "rule_dir",
        "cli_call",
        ("biological_replicates", "biological_replicates" not in samples),
        ("technical_replicates", "technical_replicates" not in samples),
        ("filter_bam_by_strand", "strandedness" not in samples),
        ("final_bam_dir", no_aligners),
        ("aligner", no_aligners),
        ("bam_sort_order", no_aligners),
        ("bam_sorter", no_aligners),
        ("deeptools_bamcoverage", no_aligners),
        ("deeptools_multibamsummary", no_aligners),
        ("deeptools_plotcorrelation", no_aligners),
        ("deeptools_qc", no_aligners),
        ("markduplicates", no_aligners),
        ("min_mapping_quality", no_aligners),
        ("only_primary_align", no_aligners),
        ("remove_blacklist", no_aligners),
        ("tx2gene_from_gtf", config.get("quantifier") != "salmon"),
        ("tximeta", config.get("quantifier") != "salmon"),
        ("deseq2", not config.get("contrasts")),
        ("deseq2_dir", not config.get("contrasts")),
        ("bigwig_dir", not config.get("create_trackhub")),
        ("trackhub_dir", not config.get("create_trackhub")),
        ("qc_dir", not config.get("create_qc_report")),
    ]
    keys = rmkeys(["samples"] + keys_to_remove, keys)
    keys = ["samples"] + keys

    ignore_values = ["", False, -1, 0, 999, "None", "none@provided.com", "yourmail@here.com", "_custom"]
    table = [(key, config[key]) for key in keys if config[key] not in ignore_values]
    table = [["\n".join(textwrap.wrap(str(cell), 112)) for cell in row] for row in table]
    table.append(("layout", tabulate([(sample, values["layout"]) for sample, values in SAMPLEDICT.items()], headers=["sample", "layout"])))

    logger.info(tabulate(table, headers=["config variable", "value"], tablefmt="pipe"))
    logger.info("\n\n")
    logger.info(f"Using the snakemake {workflow.scheduler_type} scheduler.\n")
