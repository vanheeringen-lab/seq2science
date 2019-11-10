#!/bin/bash

#generate a empty bam_file_list.txt document

for plate in $1plate*/
do
    truncate -s 0 ${plate}bam_file_list.txt
done

for plate in $1plate*/
do
    for bamfile in ${plate}cell_ID_BAMs/*.bam
    do
        echo ${bamfile} >> ${plate}bam_file_list.txt
    done
done

