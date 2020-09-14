"""
Utility functions for seq2science
"""
import os
from typing import List

import pandas as pd
import pysradb
from snakemake.io import expand


def _sample_to_idxs(df: pd.DataFrame, sample: str) -> List[int]:
    """
    Get the index/indices that belong to a run in a pysradb dataframe
    """
    if sample.startswith(("SRR", "DRR", "ERR")):
        idxs = df.index[df.run_accession == sample].tolist()
        assert len(idxs) == 1, f"sample {sample} with idxs: {idxs}"
    elif sample.startswith("SRX"):
        idxs = df.index[df.experiment_accession == sample].tolist()
        assert len(idxs) >= 1, len(idxs)
    else:
        assert False
    return idxs


def samples2metadata(samples: List[str], config: dict) -> dict:
    """
    Get the required info to continue a seq2science run from a list of samples.

    - If a sample already exists locally, we only want to know if it is paired-end or single-end.
    - If a sample does not exist locally
      - find its corresponding SRX number and all runs that belong to it,
      - check if they all have the same layout, if not, crash
      - see if we can download the runs from ena

    output:
        dict(
            "GSM1234": {"layout": "PAIRED",
                         "runs": ["SRR1234", "SRR4321"],
                         "ena_fastq_http": {
                            "SRR1234": [...],
                            "SRR4321": None
                            },
                         "ena_fastq_ftp": ["..."],

            "SRR5678": {"Layout": "SINGLE",
                        "runs": ["SRR5678"],
                        ena_fastq_http: None,
                        ena_fastq_ftp: [...],
            ...
        )
    """
    # start with empty dictionary which we fill out later
    sampledict = {sample: dict() for sample in samples}

    # fill out the sampledict for the local samples, and store the public samples for later
    public_samples = []
    for sample in samples:
        if os.path.exists(expand(f'{{fastq_dir}}/{sample}.{{fqsuffix}}.gz', **config)[0]):
            sampledict[sample]["layout"] = "SINGLE"
        elif all(os.path.exists(path) for path in expand(f'{{fastq_dir}}/{sample}_{{fqext}}.{{fqsuffix}}.gz', **config)):
            sampledict[sample]["layout"] = "PAIRED"
        elif sample.startswith(('GSM', 'SRX', 'SRR', 'ERR', 'DRR')):
            public_samples.append(sample)
        else:
            raise ValueError(f"\nsample {sample} was not found..\n"
                             f"We checked for SE file:\n"
                             f"\t{config['fastq_dir']}/{sample}.{config['fqsuffix']}.gz \n"
                             f"and for PE files:\n"
                             f"\t{config['fastq_dir']}/{sample}_{config['fqext1']}.{config['fqsuffix']}.gz \n"
                             f"\t{config['fastq_dir']}/{sample}_{config['fqext2']}.{config['fqsuffix']}.gz \n"
                             f"and since the sample did not start with either GSM, SRX, SRR, ERR, and DRR we "
                             f"couldn't find it online..\n")

    if len(public_samples) == 0:
        return sampledict

    # only continue with public samples
    db = pysradb.SRAweb()

    # all samples that are on GEO (GSM numbers), must first be converted to SRA ids (SRX numbers)
    geo_samples = [sample for sample in public_samples if sample.startswith("GSM")]

    # in sample2clean we store the (potential GEO) sample name in a SRA compliant name
    df = db.gsm_to_srx(geo_samples)
    sample2clean = dict(zip(df.experiment_alias, df.experiment_accession))

    # now add the already SRA compliant names with a reference to itself
    sample2clean.update({sample: sample for sample in public_samples if sample not in geo_samples})

    # check our samples on sra
    df = db.sra_metadata(list(sample2clean.values()), detailed=True)

    for sample, clean in sample2clean.items():
        # table indices
        idxs = _sample_to_idxs(df, clean)

        # get all runs that belong to the sample
        runs = df.loc[idxs].run_accession.tolist()
        assert len(runs) >= 1
        sampledict[sample]["runs"] = runs

        # get the layout
        layout = df.loc[idxs].library_layout.tolist()
        assert len(set(layout)) == 1, f"sample {sample} consists of mixed layouts, bad!"
        assert layout[0] in ["PAIRED", "SINGLE"], f"sample {sample} is an unclear layout, bad!"
        sampledict[sample]["layout"] = layout[0]

        # get the ena url
        sampledict[sample]["ena_fastq_http"] = dict()
        sampledict[sample]["ena_fastq_ftp"] = dict()
        for run in runs:
            if layout[0] == "SINGLE":
                sampledict[sample]["ena_fastq_http"][run] = df[df.run_accession == run].ena_fastq_http.tolist()
                sampledict[sample]["ena_fastq_ftp"][run] = df[df.run_accession == run].ena_fastq_ftp.tolist()
            elif layout[0] == "PAIRED":
                sampledict[sample]["ena_fastq_http"][run] = df[df.run_accession == run].ena_fastq_http_1.tolist() + df[
                    df.run_accession == run].ena_fastq_http_2.tolist()
                sampledict[sample]["ena_fastq_ftp"][run] = df[df.run_accession == run].ena_fastq_ftp_1.tolist() + df[
                    df.run_accession == run].ena_fastq_ftp_2.tolist()

        # if any run from a sample is not found on ENA, better be safe, and assume that sample as a whole is not on ENA
        if any(["N/A" in urls for run, urls in sampledict[sample]["ena_fastq_http"].items()]):
            sampledict[sample]["ena_fastq_http"] = None
        if any(["N/A" in urls for run, urls in sampledict[sample]["ena_fastq_ftp"].items()]):
            sampledict[sample]["ena_fastq_ftp"] = None

    return sampledict
