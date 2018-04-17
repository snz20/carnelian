#!/usr/bin/env Rscript

#' Performs 10-fold cross-validation using Random Forest classifier to check 
#' discriminatory power of functional markers selected by Carnelian. 
#'
#' @param counts_file Path of a file containing a numeric matrix of effective 
#' 					  counts of fragments assigned per functional label per sample. Note that, any 
#' 					  count matrix will do but, we advise using effective counts for more relevant 
#' 					  biological results.
#' @param sampleinfo_file Path of a tab-separated file with file with list with case-control labels of samples.
#' @param reference_group A string containing the control group name
#' @param selected_features	Path of a file containing the list of functional markers selected by Carnelian
#' @param num_folds An integer indicating number of cross-validation folds
#' @param seed	An integer for seeding the classification. Needed for reproducibility.
#' @param out_dir Path to a directory where ROC curve and classification performances will be stored. 
#' 				  Assumes the directory already exists
#' output ROC-curve and performance stats

# import required packages
library(tools)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(caret)
library(cvAUC)
library(pROC)
library(ROCR)
library(randomForest)

get_split <- function(k,n,flds){
	trainInd <- c()
	for (i in 1:n){
		if (i != k)
			trainInd <- c(trainInd,flds[[i]])
	}
	return(trainInd)
}

get_folds <- function(gp,nf){
	flds <- createFolds(gp,k=nf,list=TRUE,returnTrain=FALSE)
	return(flds)
}

draw_ROC <- function(outC, outfile){
	png(outfile,width=2200,height=2000,bg="white")
	par(mar=c(8,8,6,2)+.1,ps=36,cex=2)
	plot(outC$perf, col="grey82", lty=3, lwd=3, ann=FALSE)
	plot(outC$perf, col="darkorange", avg="vertical", lwd=10, ann=FALSE, add=TRUE)
	mtext(side = 1, text = "1 - Specificity", line = 5, cex=2)
	mtext(side = 2, text = "Sensitivity", line = 5, cex =2)
	grid(lwd=2)
	lines(x = c(0,1), y = c(0,1), col="grey82",lwd=2,lty="dotted")
	text(0.1,0.08,labels=sprintf("AUC = %.3f (95%% CI: %.2f - %.2f)", outC$cvAUC, outC$ci[1], outC$ci[2]),pos=4)
	dev.off()
}

run_classification <- function(counts,group,k,seed){
	names(counts) <- make.names(names(counts))
	accuracies <- c()
	aucs <- c()
	flds <- get_folds(group,k)
	predprobs <- list()
	predlabels <- list()
	for (i in 1:k){
		set.seed(seed)
		trainIndex <- get_split(i,k,flds)
		trainingSet <- counts[trainIndex,]
		trGrp <- group[trainIndex]
		testingSet <- counts[flds[[i]],]
		tsGrp <- group[flds[[i]]]
		modelFit <- randomForest(ClassLabel ~ .,data=trainingSet)
		
		pr <- predict(modelFit, testingSet,type='prob')
		predprobs[[i]] <- pr[,2]
		predlabels[[i]] <- testingSet$ClassLabel
		
		pred <- prediction(pr[,2],testingSet$ClassLabel)
		
		perf_AUC=performance(pred,"auc") #Calculate the AUC value
		AUC=perf_AUC@y.values[[1]]
		aucs <- c(aucs,AUC)
	}
	#print(mean(aucs))
	#print(max(aucs))
	out <- cvAUC(predprobs,predlabels,levels(group))
	outci <- ci.cvAUC(predprobs,predlabels,levels(group))
	return(c(out,outci))
}

args = commandArgs(trailingOnly=TRUE)

if (length(args)<7) {
	stop("Insufficient arguments.
	Required parameters: 
	counts_file
	sampleinfo_file
	reference_group
	selected_markers_file
	num_folds
	seed
	out_dir"
, call.=FALSE)
} else {
	counts_file <- args[1]
	countdata <- read.table(counts_file,sep='\t',header=T,row.names=1)
	countdata <- as.matrix(countdata)
	
	sampleinfo_file <- args[2]
	sampleinfo <- read.table(sampleinfo_file,sep='\t',header=T,row.names=1)
	names(sampleinfo) <- c("group","fraglen")
	
	countdata <- as.matrix(countdata[,match(rownames(sampleinfo), colnames(countdata))])
	groups <- sampleinfo$group
	ref_group <- args[3]
	groups <- relevel(groups, ref_group)
	
	selfeature_names <- as.factor(t(read.table(args[4])))
	selfeatures <- make.names(selfeature_names)
	
	dm_counts <- as.data.frame(cbind(t(countdata),groups))
	names(dm_counts)[2193] <- "ClassLabel"
	dm_counts[,"ClassLabel"] <- as.factor(groups)
	colnames(dm_counts) <- make.names(names(dm_counts))
	
	k <- as.integer(args[5])
	s <- as.integer(args[6])
	
	dmc <- dm_counts[,c(selfeatures,"ClassLabel")]
	dmc[,"ClassLabel"] <- as.factor(groups)
	out <- run_classification(dmc,groups,k,s)
	
	out_dir <- args[7]
	outpath <- file_path_as_absolute(out_dir)
	
	tablefile <- paste(outpath, 'performance.tsv', sep='/')
	x <- data.frame(flds=c(1:k,"avg"), auc=c(out$fold.AUC,out$cvAUC), ci= c(rep(NA,k), paste(toString(specify_decimal(out$ci[1],2)), toString(specify_decimal(out$ci[2],2)),sep="-")), confidence=c(rep(NA,k),out$confidence))
	write.table(x, tablefile, sep='\t', quote=FALSE, row.names=FALSE)
	
	rocfile <- paste(outpath, 'roc.png', sep='/')
	draw_ROC(out, rocfile)
}