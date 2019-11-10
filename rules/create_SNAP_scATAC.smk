rule create_SNAP_object:
    '''Create a snapobject for each BAM file, these snapobjects can be merged later using snaptools in R'''
    input:
        output_loc + '{plate}/f2q30_merged_sorted.bam'
    output:
        output_loc + '{plate}/merged_ordered.snap'
    threads: 4
    conda:
        "envs/Snaptools.yml"
    shell:
        '''
        snaptools snap-pre  \
        --input-file={input}  \
        --output-snap={output}  \
        --genome-name=hg38  \
        --genome-size="genome_files/HG38_chrom_sizes" \
        --min-mapq=30  \
        --min-flen=0  \
        --max-flen=1000  \
        --keep-chrm=TRUE  \
        --keep-single=FALSE  \
        --keep-secondary=FALSE  \
        --overwrite=True  \
        --min-cov=100  \
        --verbose=True
        '''