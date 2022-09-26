"""
- seqIO.parse(genomes) feeds contig sequences to get_masked_regions().
- get_masked_regions() returns masked regions per contif in BED3 format.
- mp.imap_unordered) runs this in parallel.
- as soon as a contig is scanned, the output is written to file.

If a contig takes longer to scan than the timeout, the scripts terminates with an error.
This prevent the code from hanging indefinitely.

Note: executing multiprocessing code inside a snakemake rule with the `run` directive
is extremely error-prone (the pool often cannot kill all its children).
"""
from multiprocessing import Pool
from Bio import SeqIO  # noqa


genome = snakemake.input.genome[0]  # noqa
mask_unsorted = snakemake.output.mask_unsorted[0]  # noqa
threads = snakemake.threads  # noqa


def get_masked_regions(contig):
    """
    Return all softmasked regions as a BED3 format string.
    """
    masked_regions = ""

    seq = contig.seq
    contig_end = len(seq)-1
    masked = False
    for i, n in enumerate(seq):
        # mark the starting position of a softmasked region
        if not masked and n.islower():
            masked = True
            start = i

        # mark end position of softmasked region (can be end of contig)
        if masked and (n.isupper() or i == contig_end):
            masked = False
            end = i

            # store softmasked region in bed3 (chr, start, end) format
            masked_regions += f"{contig.id}\t{start}\t{end}\n"  # noqa: start exists

    return masked_regions


with open(genome, "r") as genome_handle, open(mask_unsorted, "w+") as bed_handle:
    p = Pool(threads)
    contigs = p.imap_unordered(get_masked_regions, SeqIO.parse(genome_handle, "fasta"))
    p.close()
    while True:
        try:
            softmasked_regions_per_contig = contigs.next(timeout=120)
            bed_handle.write(softmasked_regions_per_contig)
        except StopIteration:
            break
