#!/bin/bash
multiBamSummary BED-file --BED cut.bed --bamfiles *.bam -o results.npz --smartLabels -p 30 --outRawCounts cut_counts.xls --scalingFactors scale.txt

while IFS="	" read -ra ADDR
do
line=${ADDR[0]}
scale=${ADDR[1]}

if [ $line = "sample" ];then

continue

fi

echo $line

bamCoverage -v -p 60 -b ${line}".bam" -o ${line}"_Crick_norm.bw" --binSize 1  --Offset 1 --samFlagInclude 128 --filterRNAstrand forward --scaleFactor $scale
bamCoverage -v -p 60 -b ${line}".bam" -o ${line}"_Watson_norm.bw" --binSize 1  --Offset 1 --samFlagInclude 128 --filterRNAstrand reverse --scaleFactor $scale


done<scale.txt
