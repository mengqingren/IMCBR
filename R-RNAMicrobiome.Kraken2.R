#!usr/bin/Rscript
###
#Author: Meng
#Usage:
#1 Public accession : Rscript R-RNAMicrobiome.Kraken2.R -a SRR14148566 -L SE -x /usr/local/data/index/hisat2_index/hg38/hg38_no_alt -D /data/mengqr/Database/Kraken2.Custom/
#2 Paired-end Fastq : Rscript R-RNAMicrobiome.Kraken2.R -f SRR14148566_1.fastq -r SRR14148566_2.fastq -P SRR14148566 -x /usr/local/data/index/hisat2_index/hg38/hg38_no_alt -D /data/mengqr/Database/Kraken2.Custom/ 
#3 Single-end Fastq : Rscript R-RNAMicrobiome.Kraken2.R -U SRR14148566.fastq -P SRR14148566 -x /usr/local/data/index/hisat2_index/hg38/hg38_no_alt -D /data/mengqr/Database/Kraken2.Custom/ 
###

#suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(plyr)))
suppressMessages(suppressWarnings(library(magrittr)))
#suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(stringr)))


option_list = list(
	make_option(c("-a", "--accession"), type="character", default=NULL,
		help="public accession", metavar="character"),
	make_option(c("-L", "--Layout"), type="character", default="PE",
		help="accession library layout PE or SE [default %default]", metavar="character"),
	make_option(c("-x", "--index"), type="character", default = NULL,
		help="the human index corresponding methods", metavar="character"),
	make_option(c("-D", "--Kraken2DB"), type="character", default = NULL,
		help="Kraken2 db path",metavar="character"),
	make_option(c("-P", "--Prefix"), type="character", default = NULL,
		help="filename",metavar="character"),
	make_option(c("-f", "--pe1"), type="character", default = NULL,
		help="forward pair end", metavar="character"),
	make_option(c("-r", "--pe2"), type="character", default = NULL,
		help="reverse pair end",metavar="character"),
	make_option(c("-U", "--se"), type="character", default = NULL,
		help="single end",metavar="character"),
	make_option(c("-s", "--savePath"), type="character", default='./',
		help="to save output dir [default %default]", metavar="character"),
	make_option(c("-t", "--threads"), type="double", default=12, 
		help="the threads used for mapping and samtools [default %default]",metavar="number"),
	make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
		help="The pipeline is used to identify microbiome from bulk RNA-seq data")#,
#	make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this help message and exit")
)

opt_parser = OptionParser(option_list=option_list)
opt <- parse_args(OptionParser(option_list=option_list))

if (any(!is.null(opt$accession),!is.null(opt$pe1),!is.null(opt$se))==F) {
	print_help(opt_parser)
	stop("Need input correct accession of fastq", call.=FALSE)
}

for (tool in c("prefetch","fastq-dump","fastp","hisat2","samtools","bedtools","kraken2")){
	text <- system(paste("which",tool),intern=T)
	if(length(text) == 0){
		print(paste("Please check the environment path of", tool, sep=" "))
	}
}

if(is.null(opt$index)){
	stop("Need input correct path of human index", call.=FALSE)
}

if(is.null(opt$Kraken2DB)){
        stop("Need input correct path of kraken2 db", call.=FALSE)
}

if(! is.null(opt$accession)){
	print("You have an accession, don not need fastq")
	prefetch_fun=paste("prefetch",opt$accession,"-O",opt$savePath)
	print(prefetch_fun)
	system(prefetch_fun)
	fastqdump_fun=paste("fastq-dump",file.path(paste(opt$savePath,opt$accession,sep="/"),paste(opt$accession,".sra",sep="")),
		"--split-3","-O",opt$savePath)
	print(fastqdump_fun)
	system(fastqdump_fun)
	if(opt$Layout == "PE"){
		fastp_fun=paste("fastp",
			"-w",opt$threads,
			"-i",file.path(opt$savePath,paste(opt$accession,"_1.fastq",sep="")),
			"-I",file.path(opt$savePath,paste(opt$accession,"_2.fastq",sep="")),
			"-o",file.path(opt$savePath,paste(opt$accession,"_1.fq",sep="")),
			"-O",file.path(opt$savePath,paste(opt$accession,"_2.fq",sep="")))
		rm_fun = paste("rm -rf",
			file.path(opt$savePath,paste(opt$accession,"_1.fastq",sep="")),
			file.path(opt$savePath,paste(opt$accession,"_2.fastq",sep="")))
		map_fun = paste("hisat2","-p",opt$threads,
			"-x",opt$index,
			"-1",file.path(opt$savePath,paste(opt$accession,"_1.fq",sep="")),
			"-2",file.path(opt$savePath,paste(opt$accession,"_2.fq",sep="")),
			"| samtools view -@",opt$threads,
			"-bS - | samtools sort - -T ./ -@",opt$threads,
			"-o",file.path(opt$savePath,paste(opt$accession,".sorted.bam",sep="")))
	}else{
		fastp_fun=paste("fastp",
			"-w",opt$threads,
			"-i",file.path(opt$savePath,paste(opt$accession,".fastq",sep="")),
			"-o",file.path(opt$savePath,paste(opt$accession,".fq",sep="")))
		#rm_fun = paste("rm -rf",file.path(opt$savePath,paste(opt$accession,".fastq",sep="")))

		map_fun = paste("hisat2","-p",opt$threads,
			"-x",opt$index,
			"-U",file.path(opt$savePath,paste(opt$accession,".fq",sep="")),
			"| samtools view -@",opt$threads,
			"-bS - | samtools sort - -T ./ -@",opt$threads,
			"-o",file.path(opt$savePath,paste(opt$accession,".sorted.bam",sep="")))
	}
	sam_fun = paste("samtools view -@",opt$threads,"-b -f 4",
                        file.path(opt$savePath,paste(opt$accession,".sorted.bam",sep="")),">",
                        file.path(opt$savePath,paste(opt$accession,".unmmaped.bam",sep="")))
        bed_fun=paste("bamToFastq -i",file.path(opt$savePath,paste(opt$accession,".unmmaped.bam",sep="")),
                        "-fq",file.path(opt$savePath,paste(opt$accession,".unmmaped.fastq",sep="")))
        kraken2_fun=paste("kraken2 --db",opt$Kraken2DB,
                        "--threads",opt$threads,
                        "--report",file.path(opt$savePath,paste(opt$accession,".report",sep="")),
                        file.path(opt$savePath,paste(opt$accession,".unmmaped.fastq",sep="")),
                        "--use-names --use-mpa-style >",
                        file.path(opt$savePath,paste(opt$accession,".out",sep="")))
#print(1)
	print(fastp_fun);system(fastp_fun)
	print(map_fun);system(map_fun)
	print(sam_fun);system(sam_fun)
	print(bed_fun);system(bed_fun)
	print(kraken2_fun);system(kraken2_fun)
}else if(all(!is.null(opt$pe1), !is.null(opt$pe1))){
	fastp_fun=paste("fastp",
			"-w",opt$threads,
			"-i",opt$pe1,
			"-I",opt$pe2,
			"-o",file.path(opt$savePath,paste(opt$Prefix,"_1.fq",sep="")),
			"-O",file.path(opt$savePath,paste(opt$Prefix,"_2.fq",sep="")))
	map_fun = paste("hisat2","-p",opt$threads,
			"-x",opt$index,
			"-1",file.path(opt$savePath,paste(opt$Prefix,"_1.fq",sep="")),
			"-2",file.path(opt$savePath,paste(opt$Prefix,"_2.fq",sep="")),
			"| samtools view -@",opt$threads,
			"-bS - | samtools sort - -T ./ -@",opt$threads,
			"-o",file.path(opt$savePath,paste(opt$Prefix,".sorted.bam",sep="")))
	sam_fun = paste("samtools view -@",opt$threads,"-b -f 4",
			file.path(opt$savePath,paste(opt$Prefix,".sorted.bam",sep="")),">", 
			file.path(opt$savePath,paste(opt$Prefix,".unmmaped.bam",sep="")))
	bed_fun=paste("bamToFastq -i",file.path(opt$savePath,paste(opt$Prefix,".unmmaped.bam",sep="")),
			"-fq",file.path(opt$savePath,paste(opt$Prefix,".unmmaped.fastq",sep="")))
	kraken2_fun=paste("kraken2 --db",opt$Kraken2DB,
			"--threads",opt$threads,
			"--report",file.path(opt$savePath,paste(opt$Prefix,".report",sep="")),
			file.path(opt$savePath,paste(opt$Prefix,".unmmaped.fastq",sep="")),
			"--use-names --use-mpa-style >",
			file.path(opt$savePath,paste(opt$Prefix,".out",sep="")))
	print(fastp_fun);system(fastp_fun)
        print(map_fun);system(map_fun)
        print(sam_fun);system(sam_fun)
        print(bed_fun);system(bed_fun)
        print(kraken2_fun);system(kraken2_fun)
}else if( ! is.null(opt$se)){
	fastp_fun=paste("fastp",
			"-w",opt$threads,
			"-i",opt$se,
			"-o",file.path(opt$savePath,paste(opt$Prefix,".fq",sep="")))
	#rm_fun = paste("rm -rf",file.path(opt$savePath,opt$se))
	map_fun = paste("hisat2","-p",opt$threads,"-x",opt$index,
			"-U",file.path(opt$savePath,paste(opt$Prefix,".fq",sep="")),
			"| samtools view -@",opt$threads,
			"-bS - | samtools sort - -T ./ -@",opt$threads,
			"-o",file.path(opt$savePath,paste(opt$Prefix,".sorted.bam",sep="")))
	sam_fun = paste("samtools view -@",opt$threads,"-b -f 4",
			file.path(opt$savePath,paste(opt$Prefix,".sorted.bam",sep="")),">", 
			file.path(opt$savePath,paste(opt$Prefix,".unmmaped.bam",sep="")))
	bed_fun=paste("bamToFastq -i",file.path(opt$savePath,paste(opt$Prefix,".unmmaped.bam",sep="")),
			"-fq",file.path(opt$savePath,paste(opt$Prefix,".unmmaped.fastq",sep="")))
	kraken2_fun=paste("kraken2 --db",opt$Kraken2DB,
			"--threads",opt$threads,
			"--report",file.path(opt$savePath,paste(opt$Prefix,".report",sep="")),
			file.path(opt$savePath,paste(opt$Prefix,".unmmaped.fastq",sep="")),
			"--use-names --use-mpa-style >",
			file.path(opt$savePath,paste(opt$Prefix,".out",sep="")))
	print(fastp_fun);system(fastp_fun)
        print(map_fun);system(map_fun)
        print(sam_fun);system(sam_fun)
        print(bed_fun);system(bed_fun)
        print(kraken2_fun);system(kraken2_fun)
}

#print(fastp_fun);print(rm_fun); print(map_fun); print(sam_fun); print(bed_fun); print(kraken2_fun)

