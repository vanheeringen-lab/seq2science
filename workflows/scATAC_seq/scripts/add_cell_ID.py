import pysam
from pysam.libchtslib import *
import os
import re
import pandas as pd
from glob import iglob
from scipy.io import mmwrite
from scipy.sparse import csr_matrix

path = os.getcwd()
if not os.path.exists(path + "/cell_ID_BAMs/"):
    os.makedirs(path + "/cell_ID_BAMs/")

#i = 0

for i, file in enumerate(snakemake.input):
    print(file)
    fn = file.split('-')[-1]
    cell_id = fn[:-6].split('_')[-0]
    samfile = pysam.AlignmentFile(file, "rb")
    sam_out = pysam.AlignmentFile(snakemake.output[i], "wb", template=samfile)
    # i += 1
    for read in samfile:
        oldQname = read.qname
        newQname = cell_id+':'+oldQname
        read.qname = newQname
        sam_out.write(read) 
    sam_out.close()
    samfile.close()