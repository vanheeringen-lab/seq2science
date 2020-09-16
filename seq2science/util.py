"""
Utility functions for seq2science
"""
import os
import time


def prep_filelock(lock_file, max_age=10):
    """
    create the directory for the lock_file if needed
    and remove locks older than the max_age (in seconds)
    """
    os.makedirs(os.path.dirname(lock_file), exist_ok=True)

    # sometimes two jobs start in parallel and try to delete at the same time
    try:
        # ignore locks that are older than the max_age
        if os.path.exists(lock_file) and \
                time.time() - os.stat(lock_file).st_mtime > max_age:
            os.unlink(lock_file)
    except FileNotFoundError:
         pass


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
