#!/bin/bash

runs=(run1 run2)
for run in ${runs[*]}
do
    echo ${run}
    cd ./${run}/
    samples=$(ls *_R1_*.fastq.gz | awk 'BEGIN{FS="_"}{print $1}' | sort | uniq)
    for sample in ${samples[*]}
    do
        echo ${sample}
	gunzip -c ${sample}*_R1_*fastq.gz >>../merged_fastq/"Sample_"${sample}"_R1.fastq"
	gunzip -c ${sample}*_R2_*fastq.gz >>../merged_fastq/"Sample_"${sample}"_R2.fastq"
    done
    cd ..
done
