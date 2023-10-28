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

  ``

