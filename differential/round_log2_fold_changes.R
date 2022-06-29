#!/usr/bin/env Rscript

# Script to round log2 fold-changes from differential expression analysis to
# the nearest 0.1.

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))

args <- parse_args(OptionParser(option_list = list(
  make_option(
    c("-e", "--experiment_type"),
    help = "Experiment type"
  ),
  make_option(
    c("-i", "--input_to_round"),
    help = "Input file with unrounded logs"
  ),
  make_option(
    c("-o", "--intermediate_output"),
    help = "Intermediate output file with rounded logs"
  )
)))

exp_type <- args$experiment_type
unrounded<- args$input_to_round
rounded  <- args$intermediate_output

unrounded_dt <- fread(file=unrounded)
# check that the file has at least one row
if( nrow( unrounded_dt )<1 ) {
  stop( paste(  unrounded, " is empty. Please check." ) )
}

# in proteomics, where we see the issue, column says <contrast>.foldChange instead of log2foldchange
# select columns
index <- unique(grep(pattern='log2foldchange', x=colnames(unrounded_dt), fixed = TRUE)	, grep(pattern='foldChange', x=colnames(unrounded_dt), fixed = TRUE)	 )
# replace NAs with 0s
for (j in index) set(unrounded_dt, which(is.na(unrounded_dt[[j]])) , j, 0)         
# round columns
for (j in index) set(unrounded_dt, j = j, value = round(unrounded_dt[[j]], digits = 1)  )

if (exp_type=='proteomics_differential') {
	# The web application requires that p-value fields are always
	# before their matching log2foldchange field (for the same group)
	# we have seen proteomics differential experiments not respecting this
	contrasts <-  tstrsplit( x=colnames(unrounded_dt)[grep(pattern='.p-value', x=colnames(unrounded_dt), fixed = TRUE)]  , split='.'  , fixed=TRUE )[[1]]
	
	reordering <- 1:ncol(unrounded_dt)
	for (i in contrasts){
		col_p <- grep(pattern=paste0(i,'.p-value'), x=colnames(unrounded_dt), fixed = TRUE)  
		col_l <- grep(pattern=paste0(i,'.log2foldchange'), x=colnames(unrounded_dt), fixed = TRUE )  
		if ( length(col_l) == 0 ) { col_l <- grep(pattern=paste0(i,'.foldChange'), x=colnames(unrounded_dt), fixed = TRUE)  } 
		# swap columns if necessary
		if( col_l < col_p){ reordering <- replace(reordering, c(col_l, col_p), reordering[c(col_p, col_l)]) } 
	}
	# reorder
	setcolorder(unrounded_dt, reordering )
} 
# Write out
write.table( unrounded_dt, file = rounded, row.names=FALSE, quote=FALSE, sep="\t" )

