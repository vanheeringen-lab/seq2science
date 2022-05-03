import os
import shutil
import sys
import trackhub


logfile = snakemake.log[0]
ALL_ASSEMBLIES = snakemake.params.ALL_ASSEMBLIES
hub = snakemake.params.hub
genomes_dir = snakemake.params.genomes_dir
out_dir = snakemake.output[0]
ORI_ASSEMBLIES = snakemake.params.ORI_ASSEMBLIES
ucsc_names = snakemake.params.ucsc_names
HAS_ANNOTATION = snakemake.params.HAS_ANNOTATION


def chmod_r(target_dir, permissions, file_permissions=None):
    """chmod -R with optionally separate permissions for directories and files"""
    if file_permissions is None:
        file_permissions = permissions

    # directories must be executable to walk through them
    os.chmod(target_dir, permissions)
    for root, dirs, files in os.walk(target_dir):
        for directory in dirs:
            os.chmod(os.path.join(root, directory), permissions)
        for file in files:
            os.chmod(os.path.join(root, file), file_permissions)


# redirect all messages to a logfile
sys.stdout = open(logfile, 'w')
sys.stderr = sys.stdout

# upload the hub
trackhub.upload.upload_hub(hub=hub, host="localhost", remote_dir=out_dir)

# actions not supported by the Trackhub package
for assembly in ALL_ASSEMBLIES:
    asmbly = ORI_ASSEMBLIES[assembly]  # no custom suffix, if present
    hub_type = "trackhub" if ucsc_names[assembly][0] else "assembly_hub"

    # copy the trix files
    if hub_type == "assembly_hub" and HAS_ANNOTATION[asmbly]:
        for ext in ["ix", "ixx"]:
            src = os.path.join(genomes_dir, assembly, f"annotation.{ext}")
            dst = os.path.join(out_dir, asmbly, f"annotation.{ext}")
            shutil.copy(src, dst)

    # add group scaling to the composite tracks
    trackdb_file = os.path.join(out_dir, ucsc_names[assembly][1], "trackDb.txt")
    with open(trackdb_file, "r") as tf:
        contents = tf.readlines()

    with open(trackdb_file, "w") as tf:
        for line in contents:
            if line.startswith("compositeTrack"):
                line = "autoScale group\n" + line
            tf.write(line)

# make the trackhub readable for everyone (not writable)
chmod_r(out_dir, 0o0755, 0o0644)
