#!usr/bin/Rscript
###
#Author: Meng
#Usage: Rscript R-Extract.Kraken2.Fastq.R --ReportFile SRR14148566.report --OutFile SRR14148566.out --UnmappedFastq SRR14148566.unmmaped.fastq --ReadsCounts 2 --SampleName SRR14148566 --savePath ./
#Rscript 
###

#suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(plyr)))
suppressMessages(suppressWarnings(library(magrittr)))
suppressMessages(suppressWarnings(library(stringr)))

option_list = list(
	make_option(c("-F", "--ReportFile"), type="character", default = NULL,
		help="Kraken2 output file .report",metavar="character"),
	make_option(c("-U", "--OutFile"), type="character", default = NULL,
		help="Kraken2 output file .out",metavar="character"),
	make_option(c("-Q", "--UnmappedFastq"), type="character", default = NULL,
		help="Unmapped Fastq with the path",metavar="character"),
	make_option(c("-R", "--ReadsCounts"), type="double", default=5, 
		help="The read counts cut off for run PathSeq [default %default]",metavar="number"),
	make_option(c("-S", "--SampleName"), type="character", default = "test",
		help="The filename for model output", metavar="character"),
	make_option(c("-P", "--savePath"), type="character", default='./',
		help="to save output dir [default %default]", metavar="character"),
	#make_option(c("-Px", "--Prefix"), type="character", default = NULL,
	#	help="filename",metavar="character"),
	make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
		help="Extract Kraken2 species fastq")
	#make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this help message and exit")
)
opt_parser = OptionParser(option_list=option_list)
opt <- parse_args(OptionParser(option_list=option_list))

#### confirm Kraken2 species
Species.count = opt$ReadsCount
data <- fread(opt$ReportFile) %>% filter(str_detect(V1,"s__")) %>% filter(str_detect(V1,"d__Bacteria")) %>% filter(V2 > Species.count) %>% mutate(Species = str_remove_all(V1,".*s__"))

data.out <- fread(opt$OutFile)
data.out <- data.out %>% filter(V1 == "C") %>% mutate(Species=str_remove_all(V3," \\(taxid.*")) %>% mutate(taxid=str_remove_all(V3,".*taxid ") %>% str_remove_all("\\)")) %>% dplyr::select(V2,V3,Species,taxid) %>% filter(Species %in% data$Species)
write.table(data.out %>% dplyr::select(V2),file=file.path(opt$savePath,paste(opt$SampleName,".Kraken2.Fastq.ID.txt",sep="")),sep="\t",quote=F,row.names=F,col.names=F)
write.table(data.out,file=file.path(opt$savePath,paste(opt$SampleName,".Kraken2.Species.ID.txt",sep="")),sep="\t",quote=F,row.names=F,col.names=T)
#### Extract Kraken2 species fastq
seqkit_exec = system("which seqkit")
if(seqkit_exec == 1){
	print("Please install Seqkit to extract species fastq files")
	system("conda install seqkit")
	stop("Please rerun the manuscript", call.=FALSE)
}
seqkit_src <- paste("seqkit grep -f ",file.path(opt$savePath,paste(opt$SampleName,".Kraken2.Fastq.ID.txt",sep=""))," ",
	opt$UnmappedFastq," > ",opt$SampleName,".Kraken2.fastq",sep="")
print(seqkit_src)
system(seqkit_src)
