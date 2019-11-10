#!/bin/bash

for i in $1*/cell_ID_BAMs; do
    truncate -s 0 ${i%cell_ID_BAMs}bam_file_list.txt
done

for bam in $(ls -l $1*/cell_ID_BAMs/*.bam | awk '{print $NF}'); do
    path=$(echo ${bam} | cut -f 2 -d/)
    echo ${bam} >> $1${path}/bam_file_list.txt
done

