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
- gatk4 Version:4.4.0.0
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
  - Refer to the web ![Kraken2](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown)
```
kraken2-build --standard --db $DBNAME
kraken2-build --standard --threads 24 --db $DBNAME
```

* PathSeq2 -> usr build
  - Refer to the web ![PathSeq](https://gatk.broadinstitute.org/hc/en-us/articles/360035889911--How-to-Run-the-Pathseq-pipeline)
```
## Create the FASTA dict
gatk CreateSequenceDictionary -R T2T.fa
gatk CreateSequenceDictionary -R microbe.fasta

## Build FASTA index
samtools faidx T2T.fa
samtools faidx microbe.fasta

## Build the host and microbe BWA index images -> This step will cost much time and memory
gatk --java-options "-Xmx80G" BwaMemIndexImageCreator -I T2T.fa
gatk --java-options "-Xmx80G -Djava.io.tmpdir=/gatk_tmp/" BwaMemIndexImageCreator -I microbe.fasta --tmp-dir /gatk_tmp/

## Generate the host k-mer library file
gatk PathSeqBuildKmers --reference T2T.fa -O T2T.hss

## Bulid taxonomy file
# ftp://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/
# ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
gatk PathSeqBuildReferenceTaxonomy -R microbe.fasta --refseq-catalog RefSeq-release220.catalog.gz --tax-dump taxdump.tar.gz -O microbe.db
```

* PathSeq2 in IMCBR
  - When run the **R-RNAMicrobiome.PathSeq2.R** , the requirement for PathSeq including **.dict**, **.fai**, **.img**, **.db** will be generated

### Quick start
#### Step1
- **R-RNAMicrobiome.Kraken2.R** -> This script contains the steps including quality control, mapping, extract unmapped reads, initial microbiome identification with Kraken2
- **Output** -> __.out, .report, .sorted.bam, .unmmaped.bam, .unmmaped.fastq__
```
## Example 1 : Public PE : SRR15115262
Rscript R-RNAMicrobiome.Kraken2.R -a SRR15115262 -L PE -x /usr/local/data/index/hisat2_index/hg38/hg38_no_alt -D /data/mengqr/Database/Kraken2.Custom/

## Example 2 : Public SE : SRR14148566
Rscript R-RNAMicrobiome.Kraken2.R -a SRR14148566 -L SE -x /usr/local/data/index/hisat2_index/hg38/hg38_no_alt -D /data/mengqr/Database/Kraken2.Custom/

## Example 3
Rscript R-RNAMicrobiome.Kraken2.R -f ${SAMPLENAME}_1.fastq -r ${SAMPLENAME}_2.fastq -P ${SAMPLENAME} -x /usr/local/data/index/hisat2_index/hg38/hg38_no_alt -D /data/mengqr/Database/Kraken2.Custom/

## Example 4
Rscript R-RNAMicrobiome.Kraken2.R -U ${SAMPLENAME}.fastq -P ${SAMPLENAME} -x /usr/local/data/index/hisat2_index/hg38/hg38_no_alt -D /data/mengqr/Database/Kraken2.Custom/
```
#### Step2
- **R-Extract.Kraken2.Fastq.R** -> This script screens out candidate microbiome with the __reads count (usr input)__ and __prevalence >= 3 samples (in script)__, then extract the candidate species fastq sequences
- **Output** -> __.Kraken2.fastq, .Kraken2.Fastq.ID.txt, .Kraken2.Species.ID.txt__
```
Rscript R-Extract.Kraken2.Fastq.R --ReportFile ${name}.report --OutFile ${SAMPLENAME}.out --UnmappedFastq ${SAMPLENAME}.unmmaped.fastq --ReadsCounts 10 --SampleName ${SAMPLENAME} --savePath ./
```

#### Step3
- **R-RNAMicrobiome.PathSeq2.R** -> Refine microbiome reads using PathSeq, integrate the microbiome composistion across the samples
- **Output** -> __.pathseq.bam, .pathseq.bam.sbi, .PathSeq.Mapped.txt, .pathseq.txt, Project.Kraken2.Species.ID.txt__
- **Importance output** -> __Project.Kraken2_PathSeq.Species_Count.csv, Project.Kraken2_PathSeq.Species_Count_Matrix.csv__
- ***Note***: This step will __cost much time and momery__
```
Rscript R-RNAMicrobiome.PathSeq2.R -C RefSeq-release220.catalog.gz -T taxdump.tar.gz -H GRCh38.fa -D ../Kraken2.Custom/ -P ./
```

#### Step4 (Optional)
- **R-BulkRNA-VOOM+SNM.R** -> VOOM+SNM: normaliztion, however, this step is optional for the datasets under specific circumstances
```
Rscript R-BulkRNA-VOOM+SNM.R -F ${Projects.Species.Count.csv} -R ${ReadsCountCutoff} -B ${BatchEffectParameters} -C ${Responsor} -M ${Metadata.csv} -S ${SAMPLENAME.VOOM_SNM} -P ./
```
#### Step5 (Optional)
- **R-Boruta.RF.R** -> Run Boruta to identify the candidate biomarkers, construct the random forest model, calculate the evaluation index AUC, F1 score 

```
conda install r-psych, r-mvtnorm, r-caret, r-PRROC, r-caTools, r-pROC, r-e1071, r-randomForest, r-Boruta, r-mlbench, r-ROSE, r-DMwR
```
```
## Example 1 => VOOM+SNM => Kfold RF
Rscript R-Boruta.RF.R -F ${Project.VOOMSNM.VOOM.SNM.Object.csv} -M ${Metadata.csv} -V TRUE -c ${Responsor} -m All -B Original -f Kfolds -S ${SAMPLENAME} -T 2 -K 5

## Example 2 => VOOM+SNM => Split RF
Rscript R-Boruta.RF.R -F ${Project.VOOMSNM.VOOM.SNM.Object.csv} -M ${Metadata.csv} -V TRUE -c ${Responsor} -m All -B Original -f Split -S ${SAMPLENAME} -T 5

## Example 3 => Count -> CPM -> Boruta -> Original -> Kfold
Rscript R-Boruta.RF.R -F ${Project.Species.Count.csv} -R ${ReadsCountCutoff} -C ${RelativeAbundanceCutoff} -M ${Metadata.csv} -V FALSE -c ${Responsor} -m Boruta -B Original -f Kfolds -S ${SAMPLENAME} -T 1 -K 5

## Example 4 => Count -> CPM -> All -> Split -> SMOTE
Rscript R-Boruta.RF.R -F ${Project.Species.Count.csv} -R ${ReadsCountCutoff} -C ${RelativeAbundanceCutoff} -M ${Metadata.csv} -V FALSE -c ${Responsor} -m All -B SMOTE -f Split -S ${SAMPLENAME} -T 5
```

### Step 6 Denoise (Optional)
This step focuses on denoise based on the provided list of microorganisms that may be contaminants in `Denoise_Species_list.txt` (__Not in IMCBR__), which can be excluded according to the needs of the USER.

**NOTE**
The instruction for these scripts can be obtained using `- h`
```
Rscript R-RNAMicrobiome.Kraken2.R -h
Rscript R-RNAMicrobiome.PathSeq2.R -h
```

### Examples
This dataset accession is PRJNA746129 in EBI, concerned about the pathogenesis of CRC liver metastasis. The final out file in __Examples__ dir.
```
## Download file, quality control, mapping, kraken2
for i in `seq 2 8`; do Rscript R-RNAMicrobiome.Kraken2.R -a SRR1511526${i} -L PE -x ~/hg38/hg38_no_alt -D /data/mengqr/Database/Kraken2.Custom/;done

## Extract candidate fastq
ls SRR1511526*.out | while read id; do name=$(basename $id | cut -d . -f 1); Rscript R-Extract.Kraken2.Fastq.R --ReportFile ${name}.report --OutFile ${name}.out --UnmappedFastq ${name}.unmmaped.fastq --ReadsCounts 10 --SampleName ${name} --savePath ./; done

## PathSeq
Rscript R-RNAMicrobiome.PathSeq2.R -C RefSeq-release220.catalog.gz -T taxdump.tar.gz -H GRCh38.fa -D ../Kraken2.Custom/ -P ./
```


















