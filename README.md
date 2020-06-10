# DEtail-seq
## Introduction
DEtail-seq is ....
### Tools Requires:
- [python](http://www.python.org/downloads/)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- [samtools](http://www.htslib.org/)
- [deeptTools](https://github.com/deeptools/deepTools)
### Python packages Requires:
- [pyBigWig](https://github.com/deeptools/pyBigWig)
- [numpy](https://numpy.org/)
- [scipy](https://www.scipy.org/)
- [statsmodels](https://www.statsmodels.org)
- [biopython](https://biopython.org/)
## Documentation:
First,duplicated reads which had same sequences for both forward and reverse reads were removed. And our own scripts RemoveSamReads.py were used to remove duplicated reads. Then,reads were aligned to the reference genome with Bowtie2 using --local settings. For visualization, the aligned reads files (BAM) were converted to single base bigWig file with 1 bp bins using bamCoverage from deepTools.Finally, poisson distribution is used for single base call peak, which is implemented by our own script DetailCallPeak.py.\
The detailed protocols were as follows:
### 1.Remove duplication reads
```
RemoveSamReads.py test_R1.fastq.gz test_R2.fastq.gz test
```

### 2.Alignment
```
bowtie2-build --threads 20 RefGenome.fasta RefGenome
bowtie2 --local --phred33 -p 20 -t -x RefGenome -1 test_dup_R1.fq.gz -2 test_dup_R2.fq.gz 2>test_align.info|samtools view -bS -1 |samtools sort -@ 10 -m 5G -l 9 -o test.sort.bam
```
### 3.Converting BAM to single base bigWig and split strand
```
samtools index test.sort.bam
bamCoverage -v -p 60 -b test.sort.bam -o test_rev.bw --minMappingQuality 30 --binSize 1  --Offset 1 --samFlagInclude 128 --filterRNAstrand forward
bamCoverage -v -p 60 -b test.sort.bam -o test_fwd.bw --minMappingQuality 30 --binSize 1  --Offset 1 --samFlagInclude 128 --filterRNAstrand reverse
```
### 4.Single base call peak
```
DetailCallPeak.py test_fwd.bw fwd test_fwd 0.01
DetailCallPeak.py test_rev.bw rev test_rev 0.01
```
## Citation
