# DEtail-seq
DEtail-seq provides a powerful approach to detect DSB ends in multiple species.  
This pipeline is used to analyze DEtail-seq data.
## Requirements
|[python](http://www.python.org/downloads/)|[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)|[bedtools](https://bedtools.readthedocs.io/en/latest/)|[samtools](http://www.htslib.org/)|[deeptTools](https://github.com/deeptools/deepTools)|
|---|---|---|---|---|
|[**pyBigWig**](https://github.com/deeptools/pyBigWig)|[**numpy**](https://numpy.org/)|[**scipy**](https://www.scipy.org/)|[**statsmodels**](https://www.statsmodels.org)|[**biopython**](https://biopython.org/)|
## Documentation
First, duplicated reads which had same sequences for both forward and reverse reads were removed. And our own scripts RemoveSamReads.py were used to remove duplicated reads. Then, reads were aligned to the reference genome with Bowtie2 using --local settings. For visualization and spliting strand, the aligned reads files (BAM) were converted to single base bigWig file with 1 bp bins using bamCoverage from deepTools. Finally, poisson distribution is used for Hotspots calling, which is implemented by our own script Hotspotcalling.py.  
  
**The detailed protocols were as follows:**
### 1.Removing duplication reads
```
$ RemoveSamReads.py test_R1.fastq.gz test_R2.fastq.gz test
```

### 2.Alignment
```
$ bowtie2-build --threads 20 RefGenome.fasta RefGenome
$ bowtie2 --local --phred33 -p 20 -t -x RefGenome -1 test_dup_R1.fq.gz -2 test_dup_R2.fq.gz\
2>test_align.info|samtools view -bS -1 |samtools sort -@ 10 -m 5G -l 9 -o test.sort.bam
```
### 3.Converting BAM to single base bigWig and splitting strand
--minMappingQuality 30 may sometimes be needed
```
$ samtools index test.sort.bam
$ bamCoverage -v -p 60 -b test.sort.bam -o test_Crick.bw --binSize 1\
--Offset 1 --samFlagInclude 128 --filterRNAstrand forward
$ bamCoverage -v -p 60 -b test.sort.bam -o test_Watson.bw --binSize 1\
--Offset 1 --samFlagInclude 128 --filterRNAstrand reverse
```
### 4.Hotspot calling
```
usage: Hotspotcalling.py [-h] [--bw BIGWIG] [--strand {fwd,rev}]
                         [--qvalue QVALUE] [--res REST] [--prefix PREFIX]

optional arguments:
  -h, --help          show this help message and exit
  --bw BIGWIG
  --strand {fwd,rev}
  --qvalue QVALUE
  --filter FILTER     A list of comma-delimited chromosome names, containing those chromosomes that 
                      should be excluded for computing the hotspotcalling and normalization.
  --res REST          Restriction sites bed file
  --norm              normalizing bw file
  --debug
  --prefix PREFIX

$ Hotspotcalling.py --bw test_Watson.bw --strand fwd --qvalue 0.01 --prefix test_Watson --norm
$ Hotspotcalling.py --bw test_Crick.bw --strand rev --qvalue 0.01 --prefix test_Crick --norm

```
If restriction enzymes are used, you may need to use the --res parameter to add restriction site information

test_Crick_basepeaks.bed and test_Watson_basepeaks.bed are the hotspot results in bed format.
test_Crick_basepeaks.xls and test_Watson_basepeaks.xls are the detailed hotspot results.

**If you want to process multiple samples in batches, please refer to DEtailPip.sh**

## Citing this work
Xu, W., Liu, C., Zhang, Z., Sun, C., Li, Q., Li, K., Jiang, H., Li, W., & Sun, Q. (2023). [DEtail-seq is an ultra-efficient and convenient method for meiotic DNA break profiling in multiple organisms.](https://pubmed.ncbi.nlm.nih.gov/36723795/) Science China. Life sciences, 66(6), 1392–1407.
