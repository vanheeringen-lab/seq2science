#!/bin/bash
# Purpose: To run kb-wrapper for RNA velocity on multiple scRNA-seq paired end fastq files stored in the same folder.

# Default values of arguments
INDEXPATH="/home/snabel/refs/genomes/GRCh38/kb_genome/velocity+erccrep_index/"
INDEXFILE="$INDEXPATH/GRCh38_index.idx"
T2G="$INDEXPATH/GRCh38_t2g.txt"
CDNAFILE="$INDEXPATH/GRCh38_cdna_t2c.txt"
INTRONFILE="$INDEXPATH/GRCh38_intron_t2c.txt"
BARCODEFILE="/home/snabel/scrna/1col_barcode_384.tab"

# Loop through arguments and process them
for arg in "$@"
do
    case $arg in
	-f|--fastqpath)
        FASTQPATH="$2"
        shift # Remove argument name from processing
        shift # Remove argument value from processing
        ;;
        -i|--idexfile)
        INDEXFILE="$2"
        shift 
        shift 
        ;;
        -g|--t2g)
        T2G="$2"
        shift
        shift
        ;;
        -b|--barcodefile)
        BARCODEFILE="$2"
        shift
        shift
        ;;
        -c1|--cdnafile)
        CDNAFILE="$2"
        shift
        shift
        ;;
        -c2|--intronfile)
        INTRONFILE="$2"
        shift
        shift
        ;;
    esac
done

echo "Running kb-wrapper for scRNA-seq with velocity, with the following settings:"
echo "# Fastq directory: $FASTQPATH"
echo "# Index path: $INDEXFILE"

if [ -z "$FASTQPATH" ];
	then
	    echo "Oops, nothing happened! This function needs the location of the .fastq files as input! Use -f or --fastqpath to specify the location of fastqs."

	else
	    for filename in ${FASTQPATH}*R1.fastq.gz; 
		do 
			id="${filename##*/}";
			id="${id/_R1.fastq.gz/}";
			output_name="${id}_output/";
			r1="$filename";
			r2="${filename/_R1.fastq.gz/}_R2.fastq.gz";
			# check if the ouput files for this plate are present.
			if [ ! -e "${output_name}counts_unfiltered/unspliced.mtx" ]; 
			then
			    echo "Running for: ${id}";
			    nice -n 10 kb count -i ${INDEXFILE} \
			    -g ${T2G} -x 0,8,16:0,0,8:1,0,0 -w ${BARCODEFILE} \
			    --overwrite --verbose --lamanno -t 40 \
			    -o ${output_name} -c1 ${CDNAFILE} -c2 ${INTRONFILE} \
			    ${r1} ${r2} >> log.out 2>&1;
			    echo "Done. Check ${output_name}";
			else 
			    echo "Checked ${output_name}. Done.";
			fi;
		done
fi
