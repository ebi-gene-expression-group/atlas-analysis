#!/usr/bin/env Rscript

# Script to create and save a SimpleList object containing one or more
# Bioconductor objects storing expression data from an Expression Atlas
# experiment. These object are SummarizedExperiment (RNA-seq data),
# ExpressionSet (1-colour microarray data), or MAList (2-colour microarray
# data).

# Get commandline arguments.
args <- commandArgs( TRUE )

# Stop if we don't have exactly two arguments.
if( length( args ) == 0 || length( args ) > 2 ) {
	# Print a usage message.
	stop( "\nUsage:\n\tcreateAtlasExperimentSummary.R <experiment accession> [ <Atlas experiment directory> ]\n\n" )

# If there's one argument, this should be the experiment accession.
} else if( length( args ) == 1 ) {
	
	# Get the accession from the first argument.
	experimentAccession <- args[ 1 ]

	# Use the current working directory as the experiment directory.
	atlasExperimentDirectory <- getwd()

# If there are two args, the first one should be the experiment accession and
# the second one the atlas experiment directory.
} else if( length( args ) == 2 ) {
	
	# Get the accession from the first argument.
	experimentAccession <- args[ 1 ]

	# Get the experiment directory from the second argument.
	atlasExperimentDirectory <- args[ 2 ]
}

# Log where we're reading and writing.
cat( paste("\nReading and writing files in", atlasExperimentDirectory, "\n\n" ) )

# Load atlasExperimentSummary package
library( atlasExperimentSummary )

# Create experiment summary SimpleList.
experimentSummary <- summarizeAtlasExperiment( 
	experimentAccession,
	atlasExperimentDirectory
)

# Create file name to save to.
experimentSummaryFile <- paste( experimentAccession, "-atlasExperimentSummary.Rdata", sep="" )
experimentSummaryFile <- file.path( atlasExperimentDirectory, experimentSummaryFile )

# Save the object to a file.
save( experimentSummary, file = experimentSummaryFile )
