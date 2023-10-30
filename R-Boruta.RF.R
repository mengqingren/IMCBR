#!/usr/bin/Rscript
#!usr/bin/Rscript
###
#Author: Meng
#Usage:
# 1 => VOOM+SNM => Kfold RF 
# Rscript ../R-Boruta.RF.R -F TumorMicrobiome.BreastCancer.PRJNA688066.HER2_ER_PR_GRADE.NAC.VOOMSNM.VOOM.SNM.Object.csv -M TumorMicrobiome.BreastCancer.PRJNA688066.txt -V TRUE -c response_to_nac -m All -B Original -f Kfolds -S test1 -T 2 -K 5
# 2 => VOOM+SNM => Split RF
# Rscript ../R-Boruta.RF.R -F TumorMicrobiome.BreastCancer.PRJNA688066.HER2_ER_PR_GRADE.NAC.VOOMSNM.VOOM.SNM.Object.csv -M TumorMicrobiome.BreastCancer.PRJNA688066.txt -V TRUE -c response_to_nac -m All -B Original -f Split -S test1 -T 5
# 3 => Count -> CPM -> Boruta -> Original -> Kfold
# Rscript ../R-Boruta.RF.R -F Treat.Tissue.RNA.BreastCancer.PRJNA688066.Species.csv -R 2 -C 1e-05 -M TumorMicrobiome.BreastCancer.PRJNA688066.txt -V FALSE -c response_to_nac -m Boruta -B Original -f Kfolds -S test1 -T 1 -K 5
# 4 => Count -> CPM -> All -> Split -> SMOTE
# Rscript ../R-Boruta.RF.R -F Treat.Tissue.RNA.BreastCancer.PRJNA688066.Species.csv -R 2 -C 1e-05 -M TumorMicrobiome.BreastCancer.PRJNA688066.txt -V FALSE -c response_to_nac -m All -B SMOTE -f Split -S test1 -T 5
#Rscript
###

suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(optparse)))

option_list = list(
	make_option(c("-F", "--CountFile"), type="character", default = NULL,
		help="Species counts file, CSV file sep with comma, column -> species, row -> samples",metavar="character"),
	make_option(c("-M", "--Metadata"), type="character", default = NULL,
		help="Corresponding metadata containing matched samples id (First column) and group information (--ColumnProject), csv file", metavar="character"),
	make_option(c("-V", "--VOOMSNM"), type="character", default = TRUE,
		help="The input count table attribute before or after VOOMSNM, FALSE -> Before, TRUE -> After", metavar="character"),
	make_option(c("-R", "--ReadsCounts"), type="double", default=5,
		help="The read counts cut off for computing abundance, only need when VOOMSNM=FALSE [default %default]",metavar="number"),
	make_option(c("-C", "--AbundanceCutoff"), type="double", default=1e-05,
		help="The abundance CPM cut off for computing species, only need when VOOMSNM=FALSE [default %default]",metavar="number"),
	make_option(c("-c", "--ColumnProject"), type="character", default = "Project",
		help="The responsor to run RF model, Label [default %default]", metavar="character"),
	make_option(c("-m", "--RFmethod"), type="character", default = "All",
		help="The species selection to run RF model, All species or Boruta selectes species [default %default]", metavar="character"),
	make_option(c("-B", "--Balance"), type="character", default = "Original",
		help="Train datasets balance method, Original -> NoBalance, SMOTE -> Balance [default %default]", metavar="character"),
	#make_option(c("-G", "--GCol"), type="character", default = NULL,help="Group colnames", metavar="character"),
	#make_option(c("-H", "--CtrlName"), type="character", default = NULL,help="Group attribute Ctrl name",metavar="character"),
	#make_option(c("-T", "--DiseaseName"), type="character", default = NULL,help="Group attribute Disease name",metavar="character"),
	make_option(c("-f", "--CVmethod"), type="character", default = "Kfolds",
		help="The train and test method to run RF model, Split -》 70%+30% or Kfolds [default %default]", metavar="character"),
	make_option(c("-S", "--SampleName"), type="character", default = "test",
		help="The filename for model output", metavar="character"),
	make_option(c("-T", "--Times"), type="double", default=1,
		help="Times to perform RF for cross validation or split 70%+30% [default %default]",metavar="number"),
	make_option(c("-K", "--Kfolds"), type="double", default=5,
		help="Kfolds for cross validation only need in CVmethod=Kfolds [default %default]",metavar="number"),
	make_option(c("-P", "--savePath"), type="character", default='./',
		help="to save output dir [default %default]", metavar="character"),
	make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
		help="Counts -> CPM -> RF model")
	#make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this help message and exit")
)
opt_parser = OptionParser(option_list=option_list)
opt <- parse_args(OptionParser(option_list=option_list))

require(data.table)
require(psych)
require(mvtnorm)
require(caret)
require(PRROC)
require(ggplot2)
require(caTools)
require(pROC)
require(dplyr)
require(e1071)
require(randomForest)
require(ROSE)
require(DMwR)
require(tidyverse)
require(magrittr)
require(Boruta)
require(mlbench)
if (is.null(opt$CountFile) | is.null(opt$Metadata)) {
	print_help(opt_parser)
	stop("Need input correct Species count file and metadata file", call.=FALSE)
}

metadata <- read.csv(opt$Metadata,row.names=1) #
if(opt$ColumnProject %in% colnames(metadata)){
	print("You have input the correct metadata file")
}else{
	print_help(opt_parser)
	stop("Need input correct metadata file containing Run and Project attribute", call.=FALSE)
}

metadata <- metadata %>% dplyr::select(c(opt$ColumnProject))
metadata$Project = metadata[[opt$ColumnProject]] %>% str_remove_all("-") %>% str_remove_all(",") %>% str_remove_all("\\(") %>% str_remove_all("\\)")
print(metadata$Project %>% unique())
metadata <- metadata %>% dplyr::select(-c(opt$ColumnProject))

RF_M <- function(data=data,column="Project",method=c("Split","Kfolds"),balance=c("Original","SMOTE"),times=5, kfold=5){
	Res_Store = list()
	Projects = unique(data[[column]]) %>% sort()
	if(length(Projects) > 2){
		stop("Need input correct dimensional labels only for two", call.=FALSE)
	}
	Res_Store[["SampleCounts"]] = table(Projects)
	data$Label = if_else(data[[column]] == Projects[1],0,1)
	data <- data %>% dplyr::select(-c(column))
	data$Label <- factor(data$Label,levels=c(0,1))
	control = trainControl(method = "cv",number = 3)

	if(method=="Split"){
		SEED <- rep(NA,times)
		OverAll <- data.frame()
		Evaluation <- data.frame()
		AUC<-rep(NA,times)
		ROC.Data <- data.frame()
		for (i in 1:times) {
			seed = round(runif(1,1,1000000))
			set.seed(seed)
			SEED[i] = seed

			index <- createDataPartition(data$Label, p=0.7, list = F)
			train_ <- data[index, ]
			test_ <- data[-index, ]
			train_[is.na(train_)] = 0
        	train_$Label=factor(train_$Label,levels = c(0,1))
        	test_[is.na(test_)] = 0
        	test_$Label=factor(test_$Label,levels = c(0,1))

			if(balance=="SMOTE"){
				print("please note that you choose the SMOTE algorithm, it may cause the error due to the samples")

        		# calculate X Y
        		Project.2 <- (train_ %>% filter(Label == 1) %>% dim())[1]
        		Project.1 <- (train_ %>% filter(Label == 0) %>% dim())[1]
	        	if (Project.2 > Project.1) {mid=(Project.2-Project.1)/Project.1
	        		TarM=0; TarN=Project.2-Project.1
	        		X =round((TarN)*100/Project.1,0)
	        		Y =round((Project.2-TarM)*100*100/(X*Project.1),0)
	        	}else{
	        		mid=(Project.1-Project.2)/Project.2
	        		TarM=0; TarN=Project.1-Project.2
	        		 X =round((TarN)*100/Project.2,0)
	        		 Y =round((Project.1+TarM)*100*100/(X*Project.2),0)
	        	}
	        	SMOTE.Sample <- DMwR::SMOTE(Label ~ ., data = train_, perc.over = X,perc.under=Y)
	        	SMOTE.Sample$Label <- factor(SMOTE.Sample$Label,levels = c(0,1))
	        	SMOTE.Sample[is.na(SMOTE.Sample)]=0
	        	fit.rf <- train(Label ~ .,data = SMOTE.Sample,method = "rf", metric="Accuracy",trControl = control)
	        	rf.pred <- predict(fit.rf, test_)
	        	c.matrix <- confusionMatrix(rf.pred,test_$Label,mode = "everything", positive="1")
	        	rf.probs = predict(fit.rf,test_,type = "prob")
        rf.ROC = roc(response = test_$Label,predictor = rf.probs$`1`,levels = levels(test_$Label))
        ROC.Data <- data.frame(X=1-rf.ROC$specificities,Y=rf.ROC$sensitivities,
                               Times=i,Source=Projects[1],Target=Projects[2]) %>%
          rbind.data.frame(ROC.Data)
        AUC[i] <- rf.ROC$auc %>% str_remove(".* ") %>% as.numeric()
        #AUC[i] <- AUC
			OverAll <- c.matrix$overall %>% data.frame() %>% magrittr::set_colnames("Score") %>% rownames_to_column("Index") %>% mutate(Times=i) %>% rbind.data.frame(OverAll)
	        	Evaluation <- c.matrix$byClass %>% data.frame() %>% magrittr::set_colnames("Score") %>% rownames_to_column("Index") %>% mutate(Times=i) %>% rbind.data.frame(Evaluation)
	        }else if(balance=="Original"){
	        	fit.rf <- train(Label ~ .,data = train_,method = "rf", metric="Accuracy",trControl = control)
	        	rf.pred <- predict(fit.rf, test_)
	        	c.matrix <- confusionMatrix(rf.pred,test_$Label,mode = "everything", positive="1")
			rf.probs = predict(fit.rf,test_,type = "prob")
			rf.ROC = roc(response = test_$Label,predictor = rf.probs$`1`,levels = levels(test_$Label))
			ROC.Data <- data.frame(X=1-rf.ROC$specificities,Y=rf.ROC$sensitivities,Times=i,Source=Projects[1],Target=Projects[2]) %>% rbind.data.frame(ROC.Data)
			AUC[i] <- rf.ROC$auc %>% str_remove(".* ") %>% as.numeric()
			#AUC[i] <- AUC
	        	OverAll <- c.matrix$overall %>% data.frame() %>% magrittr::set_colnames("Score") %>% rownames_to_column("Index") %>% mutate(Times=i) %>% rbind.data.frame(OverAll)
	        	Evaluation <- c.matrix$byClass %>% data.frame() %>% magrittr::set_colnames("Score") %>% rownames_to_column("Index") %>% mutate(Times=i) %>% rbind.data.frame(Evaluation)
	        }else{
	        	stop("Need input correct method for balance samples", call.=FALSE)
	        }
	        Res_Store[["SEED"]] = SEED
	        Res_Store[["OverAll"]] = OverAll
	        Res_Store[["Evaluation"]] = Evaluation
		Res_Store[["AUC"]] = AUC
		Res_Store[["ROC.Data"]] = ROC.Data
		}
	}else if(method=="Kfolds") {
		SEED <- rep(NA,times)
		OverAll <- data.frame()
		Evaluation <- data.frame()
		AUC <- data.frame()
		ROC.Data <- data.frame()
		for (i in 1:times) {
			seed = round(runif(1,1,1000000))
			set.seed(seed)
			SEED[i] = seed

			folds <- createFolds(data$Label,k=kfold)
			for (k in 1:kfold){
				train_ <- data[-folds[[k]],]
				test_ <- data[folds[[k]],]

				train_[is.na(train_)] = 0
        		train_$Label=factor(train_$Label,levels = c(0,1))
        		test_[is.na(test_)] = 0
        		test_$Label=factor(test_$Label,levels = c(0,1))
        		if(balance=="SMOTE"){
					print("please note that you choose the SMOTE algorithm, it may cause the error due to the samples")

        			# calculate X Y
        			Project.2 <- (train_ %>% filter(Label == 1) %>% dim())[1]
        			Project.1 <- (train_ %>% filter(Label == 0) %>% dim())[1]
	        		if (Project.2 > Project.1) {mid=(Project.2-Project.1)/Project.1
	        			TarM=0; TarN=Project.2-Project.1
	        			X =round((TarN)*100/Project.1,0)
	        			Y =round((Project.2-TarM)*100*100/(X*Project.1),0)
	        		}else{
	        			mid=(Project.1-Project.2)/Project.2
	        			TarM=0; TarN=Project.1-Project.2
	        			X =round((TarN)*100/Project.2,0)
	        			Y =round((Project.1+TarM)*100*100/(X*Project.2),0)
	        		}
	        		SMOTE.Sample <- DMwR::SMOTE(Label ~ ., data = train_, perc.over = X,perc.under=Y)
	        		SMOTE.Sample$Label <- factor(SMOTE.Sample$Label,levels = c(0,1))
	        		SMOTE.Sample[is.na(SMOTE.Sample)]=0
	        		fit.rf <- train(Label ~ .,data = SMOTE.Sample,method = "rf", metric="Accuracy",trControl = control)
	        		rf.pred <- predict(fit.rf, test_)
	        		c.matrix <- confusionMatrix(rf.pred,test_$Label,mode = "everything", positive="1")
				rf.probs = predict(fit.rf,test_,type = "prob")
				rf.ROC = roc(response = test_$Label,predictor = rf.probs$`1`,levels = levels(test_$Label))
				AUC.1 <- rf.ROC$auc %>% str_remove(".* ") %>% as.numeric()
				AUC <- data.frame(AUC=AUC.1,Times=i,Kfold=k) %>% rbind.data.frame(AUC)
				ROC.Data <- data.frame(X=1-rf.ROC$specificities,Y=rf.ROC$sensitivities,Times=i,Kfold=k,Source=Projects[1],Target=Projects[2]) %>% rbind.data.frame(ROC.Data)
	        		OverAll <- c.matrix$overall %>% data.frame() %>% magrittr::set_colnames("Score") %>% rownames_to_column("Index") %>% mutate(Times=i,Kfold=k) %>% rbind.data.frame(OverAll)
	        		Evaluation <- c.matrix$byClass %>% data.frame() %>% magrittr::set_colnames("Score") %>% rownames_to_column("Index") %>% mutate(Times=i,Kfold=k) %>% rbind.data.frame(Evaluation)
	        	}else if(balance=="Original"){
	        		fit.rf <- train(Label ~ .,data = train_,method = "rf", metric="Accuracy",trControl = control)
	        		rf.pred <- predict(fit.rf, test_)
	        		c.matrix <- confusionMatrix(rf.pred,test_$Label,mode = "everything", positive="1")
				rf.probs = predict(fit.rf,test_,type = "prob")
				rf.ROC = roc(response = test_$Label,predictor = rf.probs$`1`,levels = levels(test_$Label))
				AUC.1 <- rf.ROC$auc %>% str_remove(".* ") %>% as.numeric()
				AUC <- data.frame(AUC=AUC.1,Times=i,Kfold=k) %>% rbind.data.frame(AUC)
				ROC.Data <- data.frame(X=1-rf.ROC$specificities,Y=rf.ROC$sensitivities,Times=i,Kfold=k,Source=Projects[1],Target=Projects[2]) %>% rbind.data.frame(ROC.Data)
	        		OverAll <- c.matrix$overall %>% data.frame() %>% magrittr::set_colnames("Score") %>% rownames_to_column("Index") %>% mutate(Times=i,Kfold=k) %>% rbind.data.frame(OverAll)
	        		Evaluation <- c.matrix$byClass %>% data.frame() %>% magrittr::set_colnames("Score") %>% rownames_to_column("Index") %>% mutate(Times=i,Kfold=k) %>% rbind.data.frame(Evaluation)
	        	}else{
	        		stop("Need input correct method for balance samples", call.=FALSE)
	        	}
			}
	        Res_Store[["SEED"]] = SEED
	        Res_Store[["OverAll"]] = OverAll
	        Res_Store[["Evaluation"]] = Evaluation
		Res_Store[["AUC"]] = AUC
                Res_Store[["ROC.Data"]] = ROC.Data
		}
	}else{stop("Need input correct methods to generate train datasets and test datasets", call.=FALSE)}
	return(Res_Store)
}

if(opt$VOOMSNM){
	CountFile <- read.csv(opt$CountFile,row.names = 1) #
	CPM2 <- CountFile %>% #dplyr::select(G.S) %>%
		rownames_to_column("Run") %>%
		merge(metadata %>% rownames_to_column("Run"),by="Run") %>%
		remove_rownames() %>%
		column_to_rownames("Run")
}else{
	require(vegan)
	#print(metadata)
	print("Start transfer")
	CountFile <- read.csv(opt$CountFile,row.names = 1) #
	CountFile[CountFile < opt$ReadsCounts] = 0
	CPM <- vegan::decostand(CountFile,method = "total",MARGIN=1)
	CPM[CPM < opt$AbundanceCutoff] = 0
	P.S <- apply(CPM, 2, sum) %>% sort(decreasing = T)
	G.S <- names(P.S[P.S > 0])
	#CPM2 <- CPM %>% dplyr::select(G.S)
	#print(rownames(metadata)); print(rownames(CPM))
	#Inter.Samples <- intersect(rownames(metadata),rownames(CPM))
	#print(Inter.Samples)
	#metadata <- metadata[Inter.Samples,] %>% data.frame(check.names=F)
	#print(metadata)
	CPM2 <- CPM %>% dplyr::select(all_of(G.S)) %>% rownames_to_column("Run") %>%
		merge(metadata %>% rownames_to_column("Run"),by="Run") %>% remove_rownames() %>%
		column_to_rownames("Run")
	print("Over -> Run Count transfer CPM")
}
#print(metadata)
Projects <- CPM2[["Project"]] %>% unique() %>% sort() #Project
print(Projects)
if(length(Projects) == 2){
	if(opt$RFmethod == "All" | opt$RFmethod != "Boruta"){
		print("Feature selected without any methods -> All Species")
		CPM.F <- CPM2
		RF_RES <- RF_M(data=CPM.F,column="Project",method=opt$CVmethod,balance=opt$Balance,times=opt$Times, kfold=opt$Kfolds)
	}else if(opt$RFmethod == "Boruta"){
		CPM.F <- CPM2
		CPM.F$Label = if_else(CPM.F[["Project"]] == Projects[1],0,1)
		CPM.F <- CPM.F %>% dplyr::select(-c("Project"))
		CPM.F$Label <- factor(CPM.F$Label,levels=c(0,1))
		boruta <- Boruta(Label ~ ., data = CPM.F, doTrace = 2, maxRuns = 1000)
		saveRDS(boruta,file=file.path(opt$savePath,paste(opt$SampleName,paste0(Projects,collapse="_"),Projects[1],Projects[2],"Boruta.rds",sep=".")))
		boruta.F <- boruta$finalDecision %>% data.frame() %>% rownames_to_column() %>%
			rename(Feature=1, Decision=2) %>%
			filter(Decision == "Confirmed") %>% mutate(Project1=Projects[1],Project2=Projects[2])
		if(nrow(boruta.F) == 0){
			stop("Can not confirm the Boruta selected features", call.=FALSE)
		}
		CPM.F <- CPM.F %>% dplyr::select(c(boruta.F$Feature)) %>% mutate(Project=CPM2[["Project"]])

		RF_RES <- RF_M(data=CPM.F,column="Project",method=opt$CVmethod,balance=opt$Balance,times=opt$Times, kfold=opt$Kfolds)
		RF_RES[["Boruta"]] <- boruta.F
	}
	saveRDS(RF_RES,file=file.path(opt$savePath,paste(opt$SampleName,paste0(Projects,collapse="_"),"RF.rds",sep=".")))
}else if(length(Projects) > 2){
	MultiProjects <- list()
	for (m in 1:(length(Projects)-1)) {
		print(Projects[m])
		for(n in (m+1):length(Projects)){
			print(Projects[n])
			CPM.F <- CPM2 %>% filter(Project %in% c(Projects[m],Projects[n])) #%>% mutate(Label = if_else(Project == Projects[m],0,1)) %>% mutate(Label=factor(Label,levels=c(0,1)))
			#MultiProjects[[paste(Projects[m],Projects[n],sep="_")]] = list()
			if(opt$RFmethod == "All" | opt$RFmethod != "Boruta"){
				print("Feature selected without any methods -> All Species")
				RF_RES <- RF_M(data=CPM.F,column="Project",method=opt$CVmethod,balance=opt$Balance,times=opt$Times, kfold=opt$Kfolds)
				MultiProjects[[paste(Projects[m],Projects[n],sep="_")]] = RF_RES
			}else if(opt$RFmethod == "Boruta"){
				CPM.F2 <- CPM.F %>% mutate(Label = if_else(Project == Projects[m],0,1)) %>% mutate(Label=factor(Label,levels=c(0,1))) %>% dplyr::select(-Project)
				boruta <- Boruta(Label ~ ., data = CPM.F2, doTrace = 2, maxRuns = 1000)
				saveRDS(boruta,file=file.path(opt$savePath,paste(opt$SampleName,paste0(Projects,collapse="_"),Projects[m],Projects[n],"Boruta.rds",sep=".")))
				boruta.F <- boruta$finalDecision %>% data.frame() %>% rownames_to_column() %>%
					rename(Feature=1, Decision=2) %>%
					filter(Decision == "Confirmed") %>% mutate(Project1=Projects[1],Project2=Projects[2])
				CPM.F2 <- CPM.F2 %>% dplyr::select(c(boruta.F$Feature)) %>% mutate(Project=CPM.F[["Project"]])

				RF_RES <- RF_M(data=CPM.F,column="Project",method=opt$CVmethod,balance=opt$Balance,times=opt$Times, kfold=opt$Kfolds)
				RF_RES[["Boruta"]] <- boruta.F
				MultiProjects[[paste(Projects[m],Projects[n],sep="_")]] = RF_RES
			}
		}
	}
	saveRDS(MultiProjects,file=file.path(opt$savePath,paste(opt$SampleName,paste0(Projects,collapse="_"),"RF.rds",sep=".")))
}else{
	stop("Only one responsor label", call.=FALSE)
}

