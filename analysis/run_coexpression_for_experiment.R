#!/usr/bin/env Rscript
# Create a gene coexpression matrix for a baseline Expression Atlas experiment
# Input: undecorated tsv file, first column gene identifier, subsequent columns

# Load BiocParallel package
suppressMessages(library(clusterSeq))

# Get the commandline arguments
args <- commandArgs( TRUE )
if( length( args ) == 2 ) {
    # Input file
    expressionsFile <- args[ 1 ]

    # Destination (if not finishing with .gz, clusterSeq will append .gz to the end)
	outputPath <- args[ 2 ]
} else {
	# Print a usage message and exit.
	stop( "\nUsage:\n\trun_coexpression_for_experiment.R <expressions file: gene identifier with quartiles> <output path .gz>\n" )
}

# read file
exp <- read.delim(expressionsFile)

# Make sure there are at least three columns (gene ID, expression
# columns) -- doesn't make sense to try coexpression with only one or two
# columns of expression data.
if( ncol( exp ) < 3 ) {
    warning( "Fewer than three columns in total. Cannot do coexpression on only one data column." )

    # Quit without exit code.
    q( save="no" )
}

exp[,2:ncol(exp)] <- sapply(exp[2:ncol(exp)], function(x) sub("(^[^,]+[,][^,]+[,])([^,]+)(,.+$)", "\\2",x))  #get the middle value for each gene/tissue
expL <- sapply(exp[,2:ncol(exp)], as.numeric) # make sure the values are numeric
rownames(expL) <- exp[,1]
expL <- log(expL) # get the natural logarithm

expL[is.na(expL)] <- 0 # turn any NAs to 0

cD <- expL[rowSums(expL) > ncol(expL),] # Filter out non-expressed genes

kClust <- kCluster(cD, matrixFile= outputPath) # run kCluster function to create coexpression matrices, set to use 10 cores. It can be changed to use less or more
