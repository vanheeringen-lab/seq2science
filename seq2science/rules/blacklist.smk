"""
rules for extending/generating genome blacklist regions
"""
import re

localrules:
    setup_blacklist,
    complement_blacklist,


# the blacklist has different output depending on what it blacklists
blacklisted_filename = []
if config.get("remove_blacklist"):
    blacklisted_filename.append("encode")
if config.get("remove_mito"):
    blacklisted_filename.append("mito")
blacklisted_filename = "_".join(blacklisted_filename)


rule setup_blacklist:
    """
    Combine the encode blacklist with mitochondrial DNA depending on config.
    """
    input:
        blacklist=rules.get_genome_blacklist.output,
        sizes=rules.get_genome_support_files.output.sizes,
    output:
        expand(
            "{genome_dir}/{{assembly}}/{{assembly}}.seq2scienceblacklist_{types}.bed",
            **{**config, **{"types": blacklisted_filename}}
        ),
    params:
        config.get("remove_blacklist"),  # helps resolve changed params
        config.get("remove_mito"),  # helps resolve changed params
    run:
        newblacklist = ""
        if config.get("remove_blacklist"):
            with open(input["blacklist"][0]) as file:
                newblacklist += file.read()

        if config.get("remove_mito"):
            with open(input["sizes"][0]) as file:
                sizesfile = file.read().strip()
                for match in re.findall("chrM.*|chrm.*|MT.*", sizesfile):
                    chrm, size = match.split("\t")
                    newblacklist += f"{chrm}\t0\t{size}\n"

        with open(output[0], "w") as f:
            f.write(newblacklist)


rule complement_blacklist:
    """
    Take the complement of the blacklist. We need this complement to tell samtools
    to only keep reads that are in the complement of the blacklist.
    """
    input:
        blacklist=rules.setup_blacklist.output,
        sizes=rules.get_genome_support_files.output.sizes,
    output:
        expand(
            "{genome_dir}/{{assembly}}/{{assembly}}.seq2scienceblacklist_complement_{types}.bed",
            **{**config, **{"types": blacklisted_filename}}
        ),
    params:
        config.get("remove_blacklist"),  # helps resolve changed params
        config.get("remove_mito"),  # helps resolve changed params
    log:
        expand("{log_dir}/complement_blacklist/{{assembly}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/complement_blacklist/{{assembly}}.benchmark.txt", **config)[0]
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        sortBed -faidx {input.sizes} -i {input.blacklist} 2>> {log} |
        complementBed -i stdin -g {input.sizes} > {output} 2>> {log}
        """
