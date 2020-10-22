#!/bin/bash
genome_index="mm10" #need change
while read line
do
if [[ $line =~ ^# ]];then
        continue
fi

echo $line

RemoveSamReads.py ${line}"_combined_R1.fastq.gz" ${line}"_combined_R2.fastq.gz" $line

mkdir $line
cd $line

bowtie2 --local --phred33 -p 15 -t -x ../$genome_index -1 ../${line}"_dup_R1.fq.gz" -2 ../${line}"_dup_R2.fq.gz" 2>${line}"_align.info"|samtools sort -@ 15 -m 5G -l 9 -O BAM -o ${line}".sort.bam"
samtools index ${line}".sort.bam"
bamCoverage -v -p 60 -b ${line}".sort.bam" -o ${line}"_Crick.bw" --binSize 1  --Offset 1 --samFlagInclude 128 --filterRNAstrand forward
bamCoverage -v -p 60 -b ${line}".sort.bam" -o ${line}"_Watson.bw" --binSize 1  --Offset 1 --samFlagInclude 128 --filterRNAstrand reverse
Hotspotcalling.py --bw ${line}"_Watson.bw" --strand fwd --prefix ${line}"_Watson" --qvalue 0.01 #Watson
Hotspotcalling.py --bw ${line}"_Crick.bw" --strand rev --prefix ${line}"_Crick" --qvalue 0.01 #Crick

cd ..

done<readme
