#!/usr/bin/bash
#BSUB -J PathSeq
#BSUB -n 4
#BSUB -o /work/bio-mengqr/PathSeq
#BSUB -e /work/bio-mengqr/PathSeq.err
#BSUB -R span[ptile=40]
#BSUB -q ser
#gatk CreateSequenceDictionary -R GRCh38.fa
#samtools faidx GRCh38.fa
#samtools faidx Klebsiella_pneumoniae.fa
#samtools faidx Bacillus_megaterium.fa
#samtools microbe.fa

#gatk BwaMemIndexImageCreator -I GRCh38.fa
#gatk BwaMemIndexImageCreator -I Klebsiella_pneumoniae.fa
#gatk BwaMemIndexImageCreator -I Bacillus_megaterium.fa
#gatk BwaMemIndexImageCreator -I microbe.fa

#gatk PathSeqBuildKmers \
#    --referencePath GRCh38.fa \
#    -O GRCh38.hss
gatk --java-options "-Xmx80G" PathSeqBuildKmers     --reference GRCh38.fa     -O GRCh38.hss
gatk PathSeqBuildReferenceTaxonomy \
    -R Klebsiella_pneumoniae.fa \
    --refseq-catalog RefSeq-release220.catalog.gz \
    --tax-dump ../../../Database/Kraken2.Custom/taxonomy/taxdump.tar.gz \
    -O Klebsiella_pneumoniae.db

gatk PathSeqBuildReferenceTaxonomy \
    -R Bacillus_megaterium.fa \
    --refseq-catalog RefSeq-release220.catalog.gz \
    --tax-dump ../../../Database/Kraken2.Custom/taxonomy/taxdump.tar.gz \
    -O Bacillus_megaterium.db

#gatk PathSeqBuildReferenceTaxonomy \
#    -R microbe.fa \
#    --refseq-catalog RefSeq-release220.catalog.gz \
#    --tax-dump ../../../Database/Kraken2.Custom/taxonomy/taxdump.tar.gz \
#    -O microbe.db

