#!/usr/bin/env Rscript
#
# Script to create a gene coexpression matrix for a baseline Expression Atlas experiment.

# Load BiocParallel package.
library(clusterSeq)
#source("kmeanCluster.R")

#############
# Functions #
#############

check_file_exists <- function( filename ) {

	if( !file.exists( filename ) ) {
		stop( paste(
			"Cannot find:",
			filename
		) )
	}
}



###############################
# Script start.


# Get the commandline arguments.
args <- commandArgs( TRUE )

# Stop if we don't have the requisite number of arguments.
# For now this is 1 -- we just want the Atlas experiment directory. We can
# build the other filenames from that.
if( length( args ) == 2 ) {
	
    # Experiment accession.
    atlasExperimentAccession <- args[ 1 ]

    # Path to the Atlas experiment directory.
	atlasExperimentDirectory <- args[ 2 ]
    
    # Name of file containing expressions matrix
    #expressionsFile <- args[ 2 ]
} else {
	# Print a usage message and exit.
	stop( "\nUsage:\n\trun_coexpression_for_experiment.R <Atlas experiment directory> <expressions filename with quartiles>\n\n" )
}

message( paste( "Experiment accession is:", atlasExperimentAccession ) )

# Create path to expressions matrix file.
expressionsFile <- file.path( atlasExperimentDirectory, paste( atlasExperimentAccession, ".tsv", sep = "" ) )

# Check the directory provided exists, die if not.
check_file_exists(expressionsFile)

message( paste( "Reading expression data from file:", expressionsFile ) )

# read file
exp <- read.delim(expressionsFile)

# Make sure there are at least four columns (gene ID, gene name, expression
# columns) -- doesn't make sense to try coexpression with only one or two
# columns of expression data.
if( ncol( exp ) < 5 ) {
    warning( "Fewer than four columns in total. Cannot do coexpression on only one data column." )
    
    # Quit without exit code.
    q( save="no" )
}

exp[,3:ncol(exp)] <- sapply(exp[3:ncol(exp)], function(x) sub("(^[^,]+[,][^,]+[,])([^,]+)(,.+$)", "\\2",x))  #get the middle value for each gene/tissue
expL <- sapply(exp[,3:ncol(exp)], as.numeric) # make sure the values are numeric
rownames(expL) <- exp[,1]
expL <- log(expL) # get the natural logarithm

                                        
expL[is.na(expL)] <- 0 # turn any NAs to 0

cD <- expL[rowSums(expL) > ncol(expL),] # Filter out non-expressed genes

# preparing the output file name
outFileName <- paste( atlasExperimentAccession, "coexpressions.tsv", sep="-")

# Add full path to experiment directory to output file.
outFileName <- file.path( atlasExperimentDirectory, outFileName )

message( paste( "Coexpression matrix will be written to gzipped file:", outFileName ) )

kClust <- kCluster(cD, matrixFile= outFileName) # run kCluster function to create coexpression matrices, set to use 10 cores. It can be changed to use less or more
