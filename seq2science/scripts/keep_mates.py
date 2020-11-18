"""
Script that, after alignment, removes the information that reads are paired.
"""
from contextlib import redirect_stdout
import pysam


with open(str(snakemake.log), "w") as f:
    with redirect_stdout(f):
        paired = 1
        proper_pair = 2
        mate_unmapped = 8
        mate_reverse = 32
        first_in_pair = 64
        second_in_pair = 128

        bam_in = pysam.AlignmentFile(snakemake.input[0], "rb")
        bam_out = pysam.AlignmentFile(snakemake.output[0], "wb", template=bam_in)
        for line in bam_in:
            line.qname = line.qname + ("\\1" if line.flag & first_in_pair else "\\2")
            line.next_reference_id = 0
            line.next_reference_start = 0
            line.flag &= ~(paired + proper_pair + mate_unmapped + mate_reverse + first_in_pair + second_in_pair)

            bam_out.write(line)
