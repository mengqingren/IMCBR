#!usr/bin/Rscript
###
#Author: Meng
#Usage: Rscript ../R-BulkRNA-VOOM+SNM.R -F Treat.Tissue.RNA.BreastCancer.PRJNA688066.Species.csv -R 2 -B progesterone_receptor_status,her2_receptor_status,estrogen_receptor_status -C response_to_nac -M TumorMicrobiome.BreastCancer.PRJNA688066.txt -S TumorMicrobiome.BreastCancer.PRJNA688066.HER2_ER_PR.NAC.VOOMSNM -P ./
#Rscript 
### Note: The metadata should not contains NA 
# Note: VOOM+SNM will produce Negative values
###
###
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(optparse)))


option_list = list(
	make_option(c("-F", "--CountFile"), type="character", default = NULL,
		help="Species counts file, CSV file sep with comma, column -> Species, row -> Samples",metavar="character"),
	make_option(c("-R", "--ReadsCounts"), type="double", default=5, 
		help="The read counts cut off for computing abundance [default %default]",metavar="number"),
	make_option(c("-B", "--BatchPatameter"), type="character", default = NULL,
		help="Corresponding batcg column separated by ',' if more than one", metavar="character"),
	make_option(c("-C", "--ColumnProject"), type="character", default = "Project",
		help="The responsor to run RF model, Label [default %default]", metavar="character"),
	make_option(c("-M", "--Metadata"), type="character", default = NULL,
		help="Corresponding metadata containing matched samples id (Run), Batch information (), group information (--ColumnProject), csv file", metavar="character"),
	make_option(c("-S", "--SampleName"), type="character", default = "test",
		help="The filename for model output", metavar="character"), # sample name
	make_option(c("-P", "--savePath"), type="character", default='./',
		help="To save output dir [default %default]", metavar="character"), # save path
	make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
		help="The microbiome table predicted from RNA to run VOOM+SNM")
	#make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this help message and exit")
)
opt_parser = OptionParser(option_list=option_list)
opt <- parse_args(OptionParser(option_list=option_list))


if (is.null(opt$CountFile) | is.null(opt$Metadata)) {
	print_help(opt_parser)
	stop("Need input correct Species count file and metadata file", call.=FALSE)
}

#metadata <- read.csv(opt$Metadata)
qcMetadata <- read.csv(opt$Metadata,row.names=1) %>% na.omit()
#print(dim(qcMetadata))
if(all(opt$ColumnProject %in% colnames(qcMetadata))){
	print("You have input the correct metadata file")
}else{
	print_help(opt_parser)
	stop("Need input correct metadata file attribute", call.=FALSE)
}

require(vegan)
CountFile <- read.csv(opt$CountFile,row.names = 1) #
CountFile[CountFile < opt$ReadsCounts] = 0

P.S <- apply(CountFile, 2, sum) %>% sort(decreasing = T)
G.S <- names(P.S[P.S > 0])

CPM <- CountFile %>% dplyr::select(G.S)
#qcMetadata <- read.csv(opt$Metadata,row.names=1)
#print(dim(CPM); print(qcMetadata))
qcData <- CPM[match(rownames(qcMetadata),rownames(CPM)),]
qcMetadata <- qcMetadata[rownames(qcData),]

if(is.null(opt$BatchPatameter)){
	stop("Need input correct metadata batch parameters", call.=FALSE)
}

parameters <- (opt$BatchPatameter %>% str_split(","))[[1]]
if(all(parameters %in% colnames(qcMetadata))){
	model.formula <- as.formula(paste("~0",opt$ColumnProject,paste0(parameters,collapse="+"),sep="+"))
	bio.formula <- as.formula(paste("~",opt$ColumnProject,sep=" "))
	adj.formula <- as.formula(paste("~",paste0(parameters,collapse="+"),sep=" "))
	print(model.formula); print(bio.formula); print(adj.formula)
}else{
	stop("Need input corresponding metadata batch parameters", call.=FALSE)
}

#vsnm <- function(){
## Load packages ##
require(limma)
require(edgeR)
require(dplyr)
require(snm)
#require(doMC)
require(tibble)
require(tidyverse)

#numCores <- detectCores()
#registerDoMC(cores=numCores)

#qcMetadata <- INPUT_METADATA # ADAPT THIS AS NEEDED
#qcData <- INPUT_COUNT_DATA # ADAPT THIS AS NEEDED

# match Run id
#qcData <- qcData[match(rownames(qcMetadata),rownames(qcData)),]
#qcMetadata <- qcMetadata[rownames(qcData),]

# Set up design matrix
#covDesignNorm <- model.matrix(~0 + disease_type_consol +
#                                  host_age + # host_age should be numeric
#                                  sex, # sex should be a factor
#                                data = qcMetadata)
#covDesignNorm <- model.matrix(~0 + Project +Batch,data = qcMetadata)
covDesignNorm <- model.matrix(model.formula,data = qcMetadata)
# Check row dimensions
dim(covDesignNorm)[1] == dim(qcData)[1]

print(colnames(covDesignNorm))
# The following corrects for column names that are incompatible with downstream processing
colnames(covDesignNorm) <- gsub('([[:punct:]])|\\s+','',colnames(covDesignNorm))
print(colnames(covDesignNorm))

# Set up counts matrix
counts <- t(qcData) # DGEList object from a table of counts (rows=features, columns=samples)

# Quantile normalize and plug into voom
dge <- DGEList(counts = counts)
vdge <<- voom(dge, design = covDesignNorm, plot = TRUE, save.plot = TRUE, normalize.method="quantile")

# List biological and normalization variables in model matrices
#bio.var <- model.matrix(~disease_type_consol,data=qcMetadata)
#bio.var <- model.matrix(~Project,data=qcMetadata)
bio.var <- model.matrix(bio.formula,data=qcMetadata)
#adj.var <- model.matrix(~host_age +sex,data=qcMetadata)
#adj.var <- model.matrix(~Batch,data=qcMetadata)
adj.var <- model.matrix(adj.formula,data=qcMetadata)

colnames(bio.var) <- gsub('([[:punct:]])|\\s+','',colnames(bio.var))
colnames(adj.var) <- gsub('([[:punct:]])|\\s+','',colnames(adj.var))
print(dim(adj.var))
print(dim(bio.var))
print(dim(t(vdge$E)))
print(dim(covDesignNorm))

snmDataObjOnly <- snm(raw.dat = vdge$E, 
                      bio.var = bio.var, 
                      adj.var = adj.var, 
                      rm.adj=TRUE,
                      verbose = TRUE,
                      diagnose = TRUE)
snmData <<- t(snmDataObjOnly$norm.dat)
saveRDS(snmDataObjOnly,file=file.path(opt$savePath,paste(opt$SampleName,".VOOM.SNM.Object.Rds",sep="")))
write.csv(snmData,file=file.path(opt$savePath,paste(opt$SampleName,".VOOM.SNM.Object.csv",sep="")))

#}
# LGG HGG => Examples
if(F){
	metadata <- read.csv("GSE184941.Batch2.LGGHGG.txt") %>% dplyr::select(Run,Batch,grade) %>%
	magrittr::set_colnames(c("Run","Batch","Project")) %>% filter(Project != "") %>%
	remove_rownames() %>% column_to_rownames("Run")

	CPM <- read.csv("Tissue.RNA.LGGHGG.GSE184941.Species.csv",row.names = 1)
	CPM[CPM < 2] = 0
	qcMetadata <- metadata
	qcData <- CPM

	require(limma)
	require(edgeR)
	require(dplyr)
	require(snm)
	#require(doMC)
	require(tibble)
	require(tidyverse)
	require(vegan)
	#numCores <- detectCores()
	#registerDoMC(cores=numCores)

	qcMetadata <- INPUT_METADATA # ADAPT THIS AS NEEDED
	qcData <- INPUT_COUNT_DATA # ADAPT THIS AS NEEDED

	# match Run id
	qcData <- qcData[match(rownames(qcMetadata),rownames(qcData)),]
	qcMetadata <- qcMetadata[rownames(qcData),]

	# Set up design matrix
	#covDesignNorm <- model.matrix(~0 + disease_type_consol +
	#                                  host_age + # host_age should be numeric
	#                                  sex, # sex should be a factor
	#                                data = qcMetadata)
	covDesignNorm <- model.matrix(~0 + Project +Batch,data = qcMetadata)

	# Check row dimensions
	dim(covDesignNorm)[1] == dim(qcData)[1]

	print(colnames(covDesignNorm))
	# The following corrects for column names that are incompatible with downstream processing
	colnames(covDesignNorm) <- gsub('([[:punct:]])|\\s+','',colnames(covDesignNorm))
	print(colnames(covDesignNorm))

	# Set up counts matrix
	counts <- t(qcData) # DGEList object from a table of counts (rows=features, columns=samples)

	# Quantile normalize and plug into voom
	dge <- DGEList(counts = counts)
	vdge <<- voom(dge, design = covDesignNorm, plot = TRUE, save.plot = TRUE, normalize.method="quantile")

	# List biological and normalization variables in model matrices
	#bio.var <- model.matrix(~disease_type_consol,data=qcMetadata)
	bio.var <- model.matrix(~Project,data=qcMetadata)
	#adj.var <- model.matrix(~host_age +sex,data=qcMetadata)
	adj.var <- model.matrix(~Batch,data=qcMetadata)

	colnames(bio.var) <- gsub('([[:punct:]])|\\s+','',colnames(bio.var))
	colnames(adj.var) <- gsub('([[:punct:]])|\\s+','',colnames(adj.var))
	print(dim(adj.var))
	print(dim(bio.var))
	print(dim(t(vdge$E)))
	print(dim(covDesignNorm))

	snmDataObjOnly <- snm(raw.dat = vdge$E, 
	                      bio.var = bio.var, 
	                      adj.var = adj.var, 
	                      rm.adj=TRUE,
	                      verbose = TRUE,
	                      diagnose = TRUE)
	snmData <<- t(snmDataObjOnly$norm.dat)
	write.csv(snmData,file="Tissue.RNA.LGGHGG.GSE184941.VOOM.SNM.Species.csv")

	snmData[snmData < 0] = 0
	snmData2 <- vegan::decostand(snmData,method = "total",MARGIN=1)
	write.csv(snmData2,file="Tissue.RNA.LGGHGG.GSE184941.VOOM.SNM.Species.CPM.csv")
}
