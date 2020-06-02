#!/bin/bash
# Using getopt

abort()
{
    echo >&2 '
***************
*** ABORTED ***
***************
'
    echo "An error occurred. Exiting..." >&2
    exit 1
}

trap 'abort' 0

set -e

###############################################
## It assumes awk, bedtools and mashmap is
## available.
## We have tested this script with
## awk 4.1.3, bedtools v2.28.0 and mashmap v2.0
## on an Ubuntu system.
###############################################

threads=1
awk="awk"
bedtools="bedtools"
mashmap="mashmap"

# Argument Parsing
print_usage_and_exit () {
    echo "Usage: $0 [-j <N> =1 default] [-b <bedtools binary path> =bedtools default] [-m <mashmap binary path> =mashmap default] -a <gtf file> -g <genome fasta> -t <txome fasta> -o <output path>"
    exit 1
}

echo "****************"
echo "*** getDecoy ***"
echo "****************"
while getopts ":a:b:o:j:h:g:t:m:" opt; do
    case $opt in
        b)
            bedtools=`realpath $OPTARG`
            echo "-b <bedtools binary> = $bedtools"
            ;;
        m)
            mashmap=`realpath $OPTARG`
            echo "-m <mashmap binary> = $mashmap"
            ;;
        a)
            gtffile=`realpath $OPTARG`
            echo "-a <Annotation GTF file> = $gtffile"
            ;;
        o)
            outfolder="$OPTARG"
            echo "-o <Output files Path> = $outfolder"
            ;;
        j)
            threads="$OPTARG"
            echo "-j <Concurrency level> = $threads"
            ;;
        g)
            genomefile=`realpath $OPTARG`
            echo "-g <Genome fasta> = $genomefile"
            ;;
        t)
            txpfile=`realpath $OPTARG`
            echo "-t <Transcriptome fasta> = $txpfile"
            ;;
        h)
            print_usage_and_exit
            ;;
        \?)
            echo "Invalid option: -$OPTARG"
            print_usage_and_exit
            ;;
        :)
            echo "Option -$OPTARG requires an argument."
            print_usage_and_exit
            ;;
    esac
done

# Required arguments
if [ -z "$gtffile" -o -z "$outfolder" -o -z "$genomefile" -o -z "$txpfile" -o -z "$mashmap" -o -z "$awk" -o -z "$bedtools" -o -z "$threads" ]
then
    echo "Error: missing required argument(s)"
    print_usage_and_exit
fi

mkdir -p $outfolder
cd $outfolder

# extracting all the exonic features to mask
echo "[1/10] Extracting exonic features from the gtf"
$awk -v OFS='\t' '{if ($3=="exon") {print $1,$4,$5}}' $gtffile > exons.bed

# masking the exonic regions from the genome
echo "[2/10] Masking the genome fasta"
$bedtools maskfasta -fi $genomefile -bed exons.bed -fo reference.masked.genome.fa

# aligning the transcriptome to the masked genome
echo "[3/10] Aligning transcriptome to genome"
$mashmap -r reference.masked.genome.fa -q $txpfile -t $threads --pi 80 -s 500

# extracting the bed files from the reported alignment
echo "[4/10] Extracting intervals from mashmap alignments"
$awk -v OFS='\t' '{print $6,$8,$9}' mashmap.out | sort -k1,1 -k2,2n - > genome_found.sorted.bed

# merging the reported intervals
echo "[5/10] Merging the intervals"
$bedtools merge -i genome_found.sorted.bed > genome_found_merged.bed

# extracting relevant sequence from the genome
echo "[6/10] Extracting sequences from the genome"
$bedtools getfasta -fi reference.masked.genome.fa -bed genome_found_merged.bed -fo genome_found.fa

# concatenating the sequence at per chromsome level to extract decoy sequences
echo "[7/10] Concatenating to get decoy sequences"
$awk '{a=$0; getline;split(a, b, ":");  r[b[1]] = r[b[1]]""$0} END { for (k in r) { print k"\n"r[k] } }' genome_found.fa > decoy.fa

# concatenating decoys to transcriptome
echo "[8/10] Making gentrome"
cat $txpfile decoy.fa > gentrome.fa

# extracting the names of the decoys
echo "[9/10] Extracting decoy sequence ids"
grep ">" decoy.fa | $awk '{print substr($1,2); }' > decoys.txt

# removing extra files
echo "[10/10] Removing temporary files"
rm exons.bed reference.masked.genome.fa mashmap.out genome_found.sorted.bed genome_found_merged.bed genome_found.fa decoy.fa reference.masked.genome.fa.fai

trap : 0
echo >&2 '
**********************************************
*** DONE Processing ...
*** You can use files `$outfolder/gentrome.fa`
*** and $outfolder/decoys.txt` with
*** `salmon index`
**********************************************
'
