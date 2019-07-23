import re
import genomepy
import pandas as pd
from snakemake.utils import validate


# read the samples file
samples = pd.read_csv(config["samples"], sep='\t')
validate(samples, schema=f"../schemas/{schema}")
samples = samples.set_index('sample')
samples.index = samples.index.map(str)


# TODO: maybe move to other locations? Cleaner?
# apply workflow specific changes
if 'condition' in samples:
    samples['condition'] = samples['condition'].str.replace(" ","")

if 'peaks' in samples:
    config['peaks'] = {k: v for d in config['peaks'] for k, v in d.items()}

if 'local_lookup' in samples:
    config["local_lookup"] = samples.local_path.to_dict()

if 'assembly' in samples:
    config['assemblies'] = set(samples['assembly'])

if 'fastq_dump' in config:
    assert config['fastq_dump'] in ['split', 'splot']

# complement the config when paths are not provided
if 'result_dir' not in config:
    config['result_dir'] = f"{os.getcwd()}"
if 'genome_dir' not in config:
    config['genome_dir'] = os.path.expanduser(genomepy.functions.config.get("genome_dir"))
if 'log_dir' not in config:
    config['log_dir']    = f"{os.getcwd()}/log"

# cut off trailing slashes
for path in ['result_dir', 'genome_dir', 'log_dir']:
    config[path] = re.split("\/$", config[path])[0]


onstart:
    # TODO: download location genomepy
    if 'genomepy' in sys.modules:
        # get the genomepy settings
        config['active_plugins'] = [p.name() for p in genomepy.plugin.get_active_plugins()]

        # disable genomepy plugins
        for plugin in ['bowtie2', 'bwa', 'gmap', 'hisat2', 'minimap2']:
            genomepy.plugin.deactivate(plugin)


def onexit(config):
    if 'genomepy' in sys.modules:
        for plugin in genomepy.plugin.plugins:
            if plugin in config['active_plugins']:
                genomepy.plugin.activate(plugin)


onerror:
    onexit(config)


onsuccess:
    onexit(config)
