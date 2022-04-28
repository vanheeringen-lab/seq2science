import subprocess as sp
import sys
import trackhub


logfile = snakemake.log[0]
all_assemblies = snakemake.params.all_assemblies
hub = snakemake.params.hub
genomes_dir = snakemake.params.genomes_dir
output = snakemake.output[0]
ori_assembly = snakemake.params.ori_assembly
get_ucsc_name = snakemake.params.get_ucsc_name
has_annotation = snakemake.params.has_annotation

# redirect all messages to a logfile
sys.stdout = open(logfile, 'w')
sys.stderr = sys.stdout

# upload the hub
trackhub.upload.upload_hub(hub=hub, host="localhost", remote_dir=output)

# actions not supported by the Trackhub package
for assembly in all_assemblies:
    asmbly = ori_assembly[assembly]  # no custom suffix, if present
    hub_type = "trackhub" if get_ucsc_name[assembly][0] else "assembly_hub"

    # copy the trix files
    if hub_type == "assembly_hub" and has_annotation[asmbly]:
        for ext in ["ix", "ixx"]:
            src = f"{genomes_dir}/{assembly}/annotation.{ext}"
            dst = f"{output}/{asmbly}/annotation.{ext}"
            sp.call(f"rsync {src} {dst}", shell=True)

    # add group scaling to the composite tracks
    trackdb_file = f"{output}/{get_ucsc_name[assembly][1]}/trackDb.txt"
    with open(trackdb_file, "r") as tf:
        contents = tf.readlines()

    with open(trackdb_file, "w") as tf:
        for line in contents:
            if line.startswith("compositeTrack"):
                line = "autoScale group\n" + line
            tf.write(line)

# make the trackhub readable for everyone (not writable)
sp.call(f"chmod -R 755 {output}", shell=True)
