import requests
import re
from typing import List

import logging


def crx2downloads(crx_id):
    """
    Get the crr (run) and the corresponding download link(s) from a crx id.
    """
    # get from the crx number the accession number and the experiment url
    url = f"https://ngdc.cncb.ac.cn/search/?dbId=gsa&q={crx_id}"
    r = requests.get(url)
    if r.status_code != 200:
        return {}

    search = re.search(f"https://ngdc.cncb.ac.cn/gsa/browse/(.+)/{crx_id}", r.text)
    if search is None:
        return {}
    crx_url = search.group(0)
    cra_id = search.group(1)

    # then get the run id and url from the experiment number page
    r = requests.get(crx_url)
    if r.status_code != 200:
        return {}

    search = re.findall(f"browse/{cra_id}/(CRR\d+)", r.text)
    if search is None:
        return {}

    # and finally find the download link(s) per run and add them to the final_res dict
    final_res = {}
    for crr_id in search:
        crr_url = f"https://ngdc.cncb.ac.cn/gsa/browse/{cra_id}/{crr_id}"

        # finally find the download links that belong to the run
        r = requests.get(crr_url)
        if r.status_code != 200:
            return []

        urls = re.findall(f"https://download[^\s]+.gz", r.text)
        # remove duplicate urls but keep order
        urls = list(dict.fromkeys(urls))
        final_res[crr_id] = urls

    return final_res


def samples2metadata_gsa(samples: List[str], logger) -> dict:
    """
    Based on a list of gsa crx numbers, this function returns the layout, runs, and download links.

    output:
        dict(
            "CRX1234": {"layout": "PAIRED",
                        "runs": ["CRR1234", "CRR4321"],
                        "gsa_fastq_http": {CRR1234: [link_1, link_2], ...,

            "SRR5678": {"layout": "SINGLE",
                        "runs": ["CRR5678"],
                        gsa_fastq_http: {CRR5678: [link_1]},
            ...
        )
    """
    failed_samples = []
    sampledict = {sample: dict() for sample in samples}
    for crx_id in samples:
        rundict = crx2downloads(crx_id)
        # if nothing returned it failed
        if not rundict:
            failed_samples.append(crx_id)

        for run, urls in rundict.items():
            # if more than 2 urls we fail
            if 1 > len(urls) >= 2:
                failed_samples.append(crx_id)
                continue

            if "runs" not in sampledict[crx_id]:
                sampledict[crx_id]["runs"] = []
                sampledict[crx_id]["gsa_fastq_http"] = {}
            sampledict[crx_id]["layout"] = "PAIRED" if len(urls) == 2 else "SINGLE"
            sampledict[crx_id]["runs"].append(run)
            sampledict[crx_id]["gsa_fastq_http"][run] = urls

    if len(failed_samples) > 0:
        logger.error(f"We had trouble querying GSA. These sample(s) failed: {failed_samples}")
        os._exit(1)  # noqa

    return sampledict


print(samples2metadata_gsa(["CRX035815"], logging))



