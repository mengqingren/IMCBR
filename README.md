![Logo](Figure.Pipeline.jpg)

<h2 align="center"> IMCBR: Infer Microbiome Composition Based on RNA sequencing data

### About

IMCBR: Infer Microbiome Composition Based on RNA sequencing data

- quaity control before and after trimming with FastQC and MultiQC (**Not in IMCBR**)
- trimming adapter sequences and removing low quality sequences with Fastp
- mapping the reads onto host reference genome with Hisat2, recommend the human **T2T genome fasta**
- extract unmapped reads for subsequent microbial pathogens identification with samtools and bedtools
- initial microbial identification with Kraken2
- refined microbial confirmation with PathSeq in GATK with optional choice for microbial reference (1 user build; 2 extract candidate microbiome reference genomes from Kraken2 library, **recommand 2**)

### Dependence
## Software
- sra-tools -> Download public datasets
- fastp
- hisat2
- samtools
- bedtools
- kraken2
- seqkit
- gatk4
- R packages
- python3

### Installation
Clone the github repository into the directory for analyzing the datasets:

  `git clone https://github.com/mengqingren/IMCBR.git`

Install the dependences:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

```
conda install R==4.2 gatk4 hisat2 samtools bedtools kraken2 sra-tools
conda install r-tidyverse r-data.table r-optparse
```
### Database

* Hisat2 Index -> Example
```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz
gunzip hs1.fa.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/genes/catLiftOffGenesV1.gtf.gz
gunzip catLiftOffGenesV1.gtf.gz
mv hs1.fa T2T.fa
mv catLiftOffGenesV1.gtf T2T.gtf
extract_exons.py T2T.gtf > T2T.exon
extract_splice_sites.py T2T.gtf > T2T.ss
hisat2-build T2T.fa --ss T2T.ss --exon T2T.exon T2T
```

* Kraken2 Database -> Example
  - Refer to the web [![Kraken2](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown)]

* PathSeq2 -> usr build
```
```

* PathSeq2 in IMCBR
  - When run the **R-RNAMicrobiome.PathSeq2.R** , the requirement for PathSeq including **.dict**,**.fai**,**.img**,**.db** will be generated

### Quick start




### Test data




















