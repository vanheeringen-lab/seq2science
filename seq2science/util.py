import os

import pysradb
from snakemake.io import expand


def sample_to_idxs(df, sample):
    if sample.startswith(("SRR", "DRR", "ERR")):
        idxs = df.index[df.run_accession == sample].tolist()
        assert len(idxs) == 1, f"sample {sample} with idxs: {idxs}"
    elif sample.startswith("SRX"):
        idxs = df.index[df.experiment_accession == sample].tolist()
        assert len(idxs) >= 1, len(idxs)
    else:
        assert False
    return idxs


def samples_metadata(samples, config):
    """
    Sequencing run codes:
    SRAXXXX
    123

    1.
    S-Sra ncbi (usa)
    E-Ebi (europe)
    D-Ddbj (japan)

    2.
    R-Run (this is a guess)

    3.
    R-Run
    X-eXperiment
    S-Sample
    P-Project
    """
    sampledict = {sample: dict() for sample in samples}

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
                             f"and since the sample did not start with either GSM, SRX, SRR, ERR, and DRR we couldn't find it online..\n")

    db = pysradb.SRAweb()

    # now check GEO
    geo_samples = [sample for sample in samples if sample.startswith("GSM")]
    df = db.gsm_to_srx(geo_samples)
    sample2clean = dict(zip(df.experiment_alias, df.experiment_accession))

    df = db.sra_metadata(list(sample2clean.values()), detailed=True)

    for sample, clean in sample2clean.items():
        # table indices
        idxs = sample_to_idxs(df, clean)

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
                sampledict[sample]["ena_fastq_http"][run] = df.loc[idxs].ena_fastq_http.tolist()
                sampledict[sample]["ena_fastq_ftp"][run] = df.loc[idxs].ena_fastq_ftp.tolist()
            elif layout[0] == "PAIRED":
                sampledict[sample]["ena_fastq_http"][run] = df.loc[idxs].ena_fastq_http_1.tolist() + df.loc[
                    idxs].ena_fastq_http_2.tolist()
                sampledict[sample]["ena_fastq_ftp"][run] = df.loc[idxs].ena_fastq_ftp_1.tolist() + df.loc[
                    idxs].ena_fastq_ftp_2.tolist()

        if any(["N/A" in urls for run, urls in sampledict[sample]["ena_fastq_http"].items()]):
            sampledict[sample]["ena_fastq_http"] = None
        if any(["N/A" in urls for run, urls in sampledict[sample]["ena_fastq_ftp"].items()]):
            sampledict[sample]["ena_fastq_ftp"] = None

    return sampledict
