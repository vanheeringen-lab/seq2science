import os


def can_convert():
    """check if we can make a conversion table at all"""
    with open(snakemake.input[0]) as gtf:
        for n, line in enumerate(gtf):
            line = line.lower()
            if "gene_id" in line and "gene_name" in line:
                return True
            if n > 100:
                return False


if not can_convert():
    with open(snakemake.output[0], "w") as out:
        out.write("The gene annotation does not contain the required attribute fields 'gene_id' and 'gene_name'.\n")
    os._exit(0)  # noqa

# loop over the gtf and store the conversion in the table
table = dict()
with open(snakemake.input[0]) as gtf:
    for line in gtf:
        try:
            attributes = line.split("\t")[8].split(";")
            id, name = None, None
            for attribute in attributes:
                attribute = attribute.strip()
                if attribute.lower().startswith("gene_id"):
                    id = attribute.split(" ")[1].strip('"')
                if attribute.lower().startswith("gene_name"):
                    name = attribute.split(" ")[1].strip('"')
            if id and name:
                table[id] = name
        except IndexError:
            # skip lines that are too short/misformatted
            continue

# save the dict
with open(snakemake.output[0], "w") as out:
    for k, v in table.items():
        out.write(f"{k}\t{v}\n")
