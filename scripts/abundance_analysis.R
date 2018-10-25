#!/usr/bin/env Rscript

#' Performs differential abundance analysis on a counts matrix for two-group (case vs controls) study. 
#'
#' @param counts_file Path of a file containing a numeric matrix of effective counts of fragments assigned 
#'					 per functional label per sample. Note that, any count matrix will do but, we advise 
#'					 using effective counts for more relevant biological results.
#' @param sampleinfo_file A tab-separated file with file with list of sample ids and corresponding case-control labels of samples.
#' @param reference_group A string containing the control group name
#' @param p-value cutoff to be used for volcano plot
#' @param log-fold-change cutoff to be used for volcano plot 
#' @param out_dir Path to a directory where analysis results will be stored. Assumes the directory already exists
#' output fold-change and pvalue estimates for functional labels, volcano-plot for significant labels

# import required packages
library(tools)
library(limma)
library(Glimma)
library(edgeR)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(ggrepel)

run_analysis <- function(tpm, groups, ref_group){
	groups <- as.factor(groups) # make sure groups are factors
	groups <- relevel(groups, ref_group) # relevel so that control is the first group
	orig_names <- levels(groups) # save original group levels
	levels(groups) <- c("control","case") # rename group levels
	
	# subsetting
	thresh <- tpm > 0.5 # remove lowly expressed functional labels
	keep <- rowSums(thresh) >= 2 # select functional labels for which at least two samples have effective counts
	counts.keep <- tpm[keep,]
	
	# perform differential expression analysis using limma, voom, and edgeR
	y <- DGEList(counts.keep)
	y <- calcNormFactors(y)
	design <- model.matrix(~0+groups)
	colnames(design) <- levels(groups)
	v <- voom(y, design, plot=FALSE)
	fit <- lmFit(v)
	cont.matrix <- makeContrasts(case_vs_control = case - control, levels = design)
	fit.cont <- contrasts.fit(fit, cont.matrix)
	fit.cont <- eBayes(fit.cont)
	
	# save results table
	limma.res <- topTable(fit.cont,coef="case_vs_control",sort.by="p",n="Inf")
	return(limma.res)
}

draw_volcano_plot <- function(results, top_results, pngfile){
	g = ggplot(data=results, aes(x=logFC, y=-log10(adj.P.Val))) +
	geom_point(aes(color=group)) +
	scale_color_manual(values = c("slateblue1","grey","springgreen","tomato")) +
	theme_bw(base_size=20) + theme(legend.position="bottom") +
	geom_text_repel(
	data = top_results,
	aes(label = name),
	size = 8,
	box.padding = unit(0.35,"lines"),
	point.padding = unit(0.3,"lines")
	)
	g + labs(colour="Significance",x="Log 2 Fold Change",y = "- Log 10 (adjusted P value)")
	ggsave(file = pngfile, dpi = 600, width = 10, height = 10, units = "in")
}

args = commandArgs(trailingOnly=TRUE)

if (length(args)<6) {
	stop("Insufficient arguments.
	Required parameters: 
	counts_file
	sampleinfo_file
	reference_group
	pval_cutoff
	lfc_cutoff
	out_dir"
, call.=FALSE)
} else {
	# read counts file
	counts_file <- args[1]
	countdata <- read.table(counts_file,sep='\t',header=T,row.names=1)
	countdata <- as.matrix(countdata)
	
	sampleinfo_file <- args[2]
	sampleinfo <- read.table(sampleinfo_file,sep='\t',header=T,row.names=1)
	names(sampleinfo) <- c("group","fraglen")
	
	countdata <- as.matrix(countdata[,match(rownames(sampleinfo), colnames(countdata))])
	groups <- sampleinfo$group
	ref_group <- args[3]
	results <- run_analysis(countdata, groups, ref_group)
	
	p_cutoff <- as.numeric(args[4])
	lfc_cutoff <- as.numeric(args[5])
	
	out_dir <- args[6]
	outpath <- file_path_as_absolute(out_dir)
	
	tablefile <- paste(outpath, 'analysis_results.tsv', sep='/')
	write.table(results, tablefile, sep='\t', quote=FALSE, row.names=TRUE)
	
	results["group"] <- "NotSignificant"
	results[which(results['adj.P.Val'] < p_cutoff & abs(results['logFC']) < lfc_cutoff ),"group"] <- "Significant"
	results[which(results['adj.P.Val'] > p_cutoff & abs(results['logFC']) > lfc_cutoff ),"group"] <- "FoldChange"
	results[which(results['adj.P.Val'] < p_cutoff & abs(results['logFC']) > lfc_cutoff ),"group"] <- "Significant&FoldChange"
	results$name <- rownames(results)
	
	# select top 10 dysregulated functional labels to be annotated on the plot
	sigfc_results <- results[which(results["group"] == "Significant&FoldChange"),]
	top_peaks <- rownames(sigfc_results[with(sigfc_results, order(logFC, adj.P.Val)),][1:5,])
	top_peaks <- union(top_peaks, rownames(sigfc_results[with(sigfc_results, order(-logFC,adj.P.Val)),][1:5,]))
	top_results <- sigfc_results[top_peaks,]
	
	pngfile <- paste(outpath, 'volcano_plot.png', sep='/')
	draw_volcano_plot(results, top_results, pngfile)
	
	print(tablefile)
	print(pngfile)
}
