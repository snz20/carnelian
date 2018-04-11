#!/usr/bin/env Rscript

library(tools)

#' TPM-like effective count generation. The R implementation is based on
#' https://gist.github.com/slowkow/c6ab0348747f86e2748b 
#'
#' Reference:
#'
#'    Lior Pachter. Models for transcript quantification from RNA-Seq.
#'    arXiv:1104.3889v2 
#'    
#'    Wagner, et al. Measurement of mRNA abundance using RNA-seq data:
#'    RPKM measure is inconsistent among samples. Theory Biosci. 24 July 2012.
#'    doi:10.1007/s12064-012-0162-3
#'    
#' @param counts A numeric matrix of raw counts of
#'  fragments assigned to each functional label.
#' @param meanProteinLength A numeric vector with average protein sequence length for each label.
#' @param meanFragmentLength A numeric vector with nominal fragment lengths.
#' @return effCounts A numeric matrix normalized by library size and protein lengths.

calc_eff_counts <- function(counts, meanProteinLength, meanFragmentLength){
	# Ensure valid arguments.
	stopifnot(length(meanProteinLength) == nrow(counts))
	stopifnot(length(meanFragmentLength) == ncol(counts))
  
	# Compute effective lengths of proteins in each library.
	effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
		meanProteinLength - meanFragmentLength[i] + 1
	}))
  
	# Exclude labels with mean protein length less than the mean fragment length.
	idx <- apply(effLen, 1, function(x) min(x) > 1)
	counts <- counts[idx,]
	effLen <- effLen[idx,]
	meanProteinLength <- meanProteinLength[idx]
 
	# Process one column at a time.
	effCounts <- do.call(cbind, lapply(1:ncol(counts), function(i) {
		rate = log(counts[,i]) - log(effLen[,i])
		denom = log(sum(exp(rate)))
		exp(rate - denom + log(1e6))
	}))

	# Copy the row and column names from the original matrix.
	colnames(effCounts) <- colnames(counts)
	rownames(effCounts) <- rownames(counts)
	return(effCounts)
}

args = commandArgs(trailingOnly=TRUE)

if (length(args)<4) {
	stop("Insufficient arguments.
	Required parameters: 
	counts_file
	gs_file
	mapping_file
	out_dir"
	, call.=FALSE)
} else {
	counts_file <- args[1]
	countdata <- read.table(counts_file,sep='\t',header=T,row.names=1)
	countdata <- as.matrix(countdata)

	gsfile <- args[2]
	gsinfo <- read.table(gsfile,sep='\t',header=T,row.names=1)
	names(gsinfo) = "genelen"	

	samplefile <- args[3]
	sampleinfo <- read.table(samplefile,sep='\t',header=T,row.names=1)
	names(sampleinfo) <- c("group","fraglen")

	countdata <- as.matrix(countdata[,match(rownames(sampleinfo), colnames(countdata))])
	countdata <- as.matrix(countdata[match(rownames(gsinfo), rownames(countdata)),])
	eff_counts <- calc_eff_counts(countdata, gsinfo$genelen, sampleinfo$fraglen)

	out_dir <- args[4]
	outpath <- file_path_as_absolute(out_dir)

	outfile <- paste(outpath,'effective_counts.tsv',sep='/')
	print(outfile)
        write.table(eff_counts,outfile,sep='\t',quote=FALSE,row.names=TRUE)
}
