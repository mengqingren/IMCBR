#!usr/bin/Rscript
###
#Author: Meng
#Usage:
#Rscript 
###

#suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(tibble)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(plyr)))
suppressMessages(suppressWarnings(library(magrittr)))
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(stringr)))
option_list = list(
	make_option(c("-C", "--Catalog"), type="character", default = "RefSeq-release220.catalog.gz",
		help="NCBI cataloh file for run PathSeq according to the GATK4 instruction, ftp://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/",metavar="character"),
	make_option(c("-T", "--Taxonomy"), type="character", default = "taxdump.tar.gz",
		help="NCBI taxonomy file for run PathSeq according to the GATK4 instruction, ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",metavar="character"),
	make_option(c("-H", "--GRCh38"), type="character", default="./GRCh38.fa", 
		help="The host genome [default %default]",metavar="character"),
	make_option(c("-D", "--Kraken2DB"), type="character", default = NULL,
		help="Kraken2 db path",metavar="character"),
	make_option(c("-P", "--savePath"), type="character", default='./',
		help="The input and output dir [default %default]", metavar="character"),
	make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
		help="Run PathSeq and generate count table")
)
opt_parser = OptionParser(option_list=option_list)
opt <- parse_args(OptionParser(option_list=option_list))

files <- list.files(path=opt$savePath,pattern=".Kraken2.Species.ID.txt")

if(length(files) == 0){
	stop("Can not find Kraken2.Species.ID.txt file,Please note the path containing the file", call.=FALSE)
}
                         
Species_Taxid = data.frame()
for (file in files) {
	Species_Taxid <- fread(file) %>% dplyr::select(Species,taxid) %>% arrange(Species) %>% unique() %>% 
	rbind.data.frame(Species_Taxid)
}
Species_Taxid <- Species_Taxid %>% group_by(Species,taxid) %>% dplyr::summarise(Count=n()) %>% filter(Count>=3) %>% dplyr::select(Species,taxid)
write.table(Species_Taxid,file="Project.Kraken2.Species.ID.txt",row.names=F,sep="\t",col.names=T,quote=F)
if(! file.exists(opt$GRCh38)){
	stop("Cannot find host reference", call.=FALSE)
}

gatk_exec = system("which gatk")
if(gatk_exec == 1){
	print("Cannot find the gatk4, Please install gatk4 and samtools in gatk4 environment to run PathSeq as following command")
	print("conda install gatk4")
	#print("conda activate GATK4")
	print("conda install samtools")
	#print("conda deactivate")
	#system("conda install gatk4")
	#system("conda activate GATK4")
	#system("conda install samtools")
	#system("conda deactivate")
	stop("Please rerun the manuscript", call.=FALSE)
}else{
	print("GATK4 has been found for PathSeq")
}

#### Extract microbiota referrence fasta
if(! file.exists("microbe.fa")){
	Extract_exec = paste("python Extract.Microbe.Fasta.py --Kraken2DB ", opt$Kraken2DB ," --TaxidFile Project.Kraken2.Species.ID.txt --Output microbe.fa --OutputRef Ref_Microbe.txt",sep="")
	print(Extract_exec)
	system(Extract_exec)
	#stop("Cannot find microbe.fa reference, plaease note the KrakenDB path", call.=FALSE)
} 
if(! file.exists("microbe.fa")){
stop("Cannot generate and find the microbe.fa", call.=FALSE)
}
print("Begin PathSeq")
#### ln host reference
print("Build Soft Linking")
if( ! file.exists("GRCh38.fa")){
	system(paste("ln -s ",opt$GRCh38," GRCh38.fa"))
}

#system("conda activate GATK4")
#### Build Reference Dict
if(! file.exists("GRCh38.dict")){
	system(paste("gatk CreateSequenceDictionary -R GRCh38.fa",sep=""))
}
if(! file.exists("microbe.dict")){
	system(paste("gatk CreateSequenceDictionary -R microbe.fa --java-options '-Xmx120G'",sep=""))
}
#### Build Reference faidx
if(! file.exists("GRCh38.fa.fai")){
	system(paste("samtools faidx GRCh38.fa",sep=""))
}
system(paste("samtools faidx microbe.fa",sep=""))

#### Build Reference bwa index
print("Please note this step will cost much time and memory")
if(! file.exists("GRCh38.fa.img")){
	system(paste("gatk BwaMemIndexImageCreator -I GRCh38.fa",sep=""))
}
#system(paste("gatk BwaMemIndexImageCreator -I microbe.fa",sep=""))
if( ! file.exists("microbe.fa.img")){
	system(paste("gatk --java-options '-Djava.io.tmpdir=./' BwaMemIndexImageCreator -I microbe.fa --tmp-dir ./",sep=""))
}
#### Build Host Reference Kmer
if(! file.exists("GRCh38.hss")){
	system(paste("gatk --java-options '-Xmx120G' PathSeqBuildKmers --reference GRCh38.fa -O GRCh38.hss",sep=""))
}
#### Build Microbial DB 
if(! file.exists("microbe.db")){
	system(paste("gatk --java-options '-Xmx150G' PathSeqBuildReferenceTaxonomy -R microbe.fa --refseq-catalog ",opt$Catalog," --tax-dump ",opt$Taxonomy," -O microbe.db",sep=""))
}
#### Run
files <- list.files(path=opt$savePath,pattern=".Kraken2.fastq")
for (file in files){
	filename = file %>% str_remove_all("\\..*")
	### fastqtobam
	system(paste("gatk FastqToSam --FASTQ ",file," --OUTPUT ",filename,".bam --SAMPLE_NAME ",filename,sep=""))
	### Run PathSeq
	PathSeq_exec = paste("gatk --java-options '-Xmx300G' PathSeqPipelineSpark --input ",filename, ".bam",
		" --filter-bwa-image GRCh38.fa.img --kmer-file GRCh38.hss --min-clipped-read-length 40 ",
		"--microbe-bwa-image microbe.fa.img ",
		"--microbe-dict  microbe.dict ",
		"--taxonomy-file microbe.db ",
		"--output ",filename,".pathseq.bam ", 
		"--scores-output ", filename, ".pathseq.txt ", 
		"--read-filter WellformedReadFilter ",
		"--divide-by-genome-length true",sep="")
	print(PathSeq_exec)
	system(PathSeq_exec)
}

files <- list.files(path=opt$savePath,pattern=".pathseq.bam$")
Species_Count <- data.frame()
for (file in files){
	filename = file %>% str_remove_all("\\..*")
	kraken2_taxid = fread(paste(filename,".Kraken2.Species.ID.txt",sep="")) %>% dplyr::select(V2,Species,taxid) %>% 
	magrittr::set_colnames(c("ID","Species","Taxid")) %>% mutate(Species=str_replace_all(Species," ","_"))
	system(paste("samtools view -F 4 ",file," | cut -f 1,3 > ", filename, ".PathSeq.Mapped.txt",sep=""))
	pathseq_taxid = fread("Ref_Microbe.txt",header=F,col.names=c("RefID","Species","Taxid")) %>%
	mutate(Species=str_remove_all(Species,"\\|.*")) %>% mutate(RefID=str_remove_all(RefID,">"))
	pathseq_map = fread(paste(filename, ".PathSeq.Mapped.txt",sep=""),header=F,col.names=c("ID","RefID")) %>% merge(pathseq_taxid,by="RefID")
	kraken2_pathseq <- merge(pathseq_map,kraken2_taxid,by=c("ID","Species","Taxid"))
	write.csv(kraken2_pathseq,file=paste(filename,".kraken2_pathseq.id_taxonomy.csv",sep=""),quote=F,row.names=F)
	Species_Count <- kraken2_pathseq %>% group_by(Species,Taxid) %>% dplyr::summarise(Count=n()) %>% mutate(Sample=filename) %>% 
	rbind.data.frame(Species_Count)
}
write.csv(Species_Count,file="Project.Kraken2_PathSeq.Species_Count.csv",row.names=F,quote=F)

Species_Count_Matrix <- Species_Count %>% dplyr::select(-Taxid) %>% remove_rownames() %>% spread(Species,Count,fill=0)
write.csv(Species_Count_Matrix,file="Project.Kraken2_PathSeq.Species_Count_Matrix.csv",row.names=F,quote=F)
