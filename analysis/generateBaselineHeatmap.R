#!/usr/bin/env Rscript
#
# Script to create a heatmap for a baseline Expression Atlas experiment.

# Load atlasConfigParser package to read XML config.
suppressMessages( library( atlasConfigParser ) )
# Load genefilter package to get rowVars function.
suppressMessages( library( genefilter ) )
# Load gplots package for heatmap.2 function.
suppressMessages( library( gplots ) )
# Load RColorBrewer package for nice colours.
suppressMessages( library( RColorBrewer ) )

# Get the commandline arguments.
args <- commandArgs( TRUE )

# Stop if we don't have the requisite number of arguments.
# For now this is 1 -- we just want the Atlas experiment directory. We can
# build the other filenames from that.
if( length( args ) == 1 ) {
	# Otherwise, collect the path to the Atlas experiment directory.
	atlasExperimentDirectory <- args[ 1 ]
} else if( length( args ) == 2 ) {
	atlasExperimentDirectory <- args[ 1 ]
	species <- args[ 2 ]
} else {
	# Print a usage message and exit.
	stop( "\nUsage:\n\tgenerateBaselineHeatmap.R <Atlas experiment directory> [<species>]\n\n" )
}

# Check the directory provided exists, die if not.
if( !file.exists( atlasExperimentDirectory ) ) {
	stop( paste( 
		"The experiment directory provided does not exist -- cannot continue.\nDirectory provided was:\n\t", 
		atlasExperimentDirectory, 
		sep = "" 
	) )
}

# Build the XML config file name and the FPKM matrix filename.
# Need to get the experiment accession out of the path passed.
experimentAccession <- basename( atlasExperimentDirectory )

# Log accession.
message( paste( "Experiment accession is:", experimentAccession ) )

# Log species if there is one.
if( exists( "species" ) ) {
	message( paste( "Species is:", species ) )
}

# XML config file.
experimentConfigFile <- file.path( atlasExperimentDirectory, paste( experimentAccession, "-configuration.xml", sep="" ) )
# Check it exists.
if( !file.exists( experimentConfigFile ) ) {
	stop( paste(
		"Cannot find XML config file for experiment",
		experimentAccession,
		"at:\n",
		experimentConfigFile
	) )
}

# If species was provided, append to accession.
experimentAccessionAndSpecies <- experimentAccession
if( exists( "species" ) ) {
	experimentAccessionAndSpecies <- paste( experimentAccessionAndSpecies, species, sep="_" )
}

# FPKMs matrix file.
fpkmsMatrixFile <- file.path( atlasExperimentDirectory, paste( experimentAccessionAndSpecies, ".tsv", sep="" ) )
# Check it exists.
if( !file.exists( fpkmsMatrixFile ) ) {
	stop( paste(
		"Cannot find FPKMs matrix file for experiment",
		experimentAccession,
		"at:\n",
		fpkmsMatrixFile
	) )
}

# Parse XML config to a list.
message( paste( "Reading experiment XML config from", experimentConfigFile, "..." ) )
experimentConfigList <- parseAtlasXML( experimentConfigFile )
message( "Successfully read experiment XML config" )

# Get the AssayGroup objects from the experiment config ready to use later.
# First get list of analytics.
experimentAnalytics <- experimentConfigList$allAnalytics
# Now get the RNA-seq analytics list (there should be the only analytics element).
rnaseqAnalytics <- experimentAnalytics$rnaseq
# Now get the list of RNA-seq assay groups.
rnaseqAssayGroups <- assay_groups( rnaseqAnalytics )

# Check that all assay groups have a label, we cannot continue without one.
message( "Verifying that all assay groups in XML config have labels..." )
invisible( lapply( rnaseqAssayGroups, function( assayGroup) {
	if( length( assay_group_label( assayGroup ) ) == 0 ) {
		stop( paste( "Assay group", assay_group_id( assayGroup ), "does not have a label. Cannot continue." ) )
	}
} ) )
message( "All assay groups have labels." )

# Read in the FPKMs.
message( paste( "Reading FPKMs from", fpkmsMatrixFile, "..." ) )
fpkmsDataFrame <- read.delim( fpkmsMatrixFile, stringsAsFactors=FALSE, header=TRUE )
message( "Successfully read FPKMs" )

# Assign gene IDs as row names.
rownames( fpkmsDataFrame ) <- fpkmsDataFrame$Gene.ID
# Remove the Gene.ID column.
fpkmsDataFrame$Gene.ID <- NULL

# Create data frame of just Ensembl IDs and gene names. We have to use Ensembl
# IDs as row names in R because it doesn't allow non-unique row names.
geneIDsToGeneNames <- data.frame( Gene.Name = fpkmsDataFrame$Gene.Name, stringsAsFactors=FALSE )
# Add the Ensembl gene IDs as the row names.
rownames( geneIDsToGeneNames ) <- rownames( fpkmsDataFrame )

# Replace any empty gene names with the gene ID.
emptyGeneNameIdxs <- which( geneIDsToGeneNames == "" )
geneIDsToGeneNames[ emptyGeneNameIdxs , ] <- rownames( geneIDsToGeneNames )[ emptyGeneNameIdxs ]

# Now remove the Gene.Name column from the data frame.
fpkmsDataFrame$Gene.Name <- NULL

# The FPKMs data frame can contain non-numeric values, such as "LOWDATA", which
# come from Cufflinks. We treat this as missing data, and in order for R to
# understand that we have to convert all non-numeric values to NA. The
# "as.numeric" function does this for us. Here we apply it to each column of
# the FPKMs data frame, to create a new one, without any non-numeric values.
message( "Converting all values to numeric..." )
fpkmsNumeric <- data.frame( lapply( fpkmsDataFrame, function( x ) {
	# Suppress warnings about coercion to NA -- that's why we're using this
	# function anyway!
	suppressWarnings( as.numeric( x ) )
} ) )
message( "Numeric conversion complete" )

# Conversion to numeric-only has removed the row names, so we have to add them back.
rownames( fpkmsNumeric ) <- rownames( fpkmsDataFrame )

# To create the heatmap, we need to take the top 100 most variable genes. To do
# this, we will use the "rowVars" function from the genefilter package, which
# calculates the variance of each row of a data frame.
rowVariances <- rowVars( fpkmsNumeric )

# Sort FPKMs by the row variances in descending order (largest first).
fpkmsNumeric <- fpkmsNumeric[ order( rowVariances, decreasing = TRUE ) , ]

# Get the FPKMs of the top 100 most variable genes.
top100geneFPKMs <- fpkmsNumeric[ 1:100 , ]

# Get the gene names for the top 100 gene IDs, to use as labels for the
# heatmape rows.
top100geneNames <- geneIDsToGeneNames[ rownames( top100geneFPKMs ) , ]

# Scale and center the expression levels using Z-score transformation, so they
# have mean 0 and standard deviation 1. This makes the clustering and heatmap
# colours look better than leaving them as they are (very long tail).
top100geneFPKMs <- t( scale( t( top100geneFPKMs ))) 

# Get the assay group labels to use as the labels for the heatmap columns.
assayGroupLabels <- sapply( colnames( top100geneFPKMs ), function( assayGroupID ) {
	assayGroup <- rnaseqAssayGroups[[ assayGroupID ]]
	assayGroupLabel <- assay_group_label( assayGroup )
} )

# Make the heatmap filename.
heatmapFilename <- paste( experimentAccession, "-heatmap.pdf", sep="" )
# Prepend path to experiment directory.
heatmapFilename <- file.path( atlasExperimentDirectory, heatmapFilename )

# Some nice colours.
colours <- colorRampPalette( brewer.pal( 9, "Blues" ) )( 100 )

# Make the heatmap width a function of the number of assay groups. If the
# columns are too narrow it's impossible to read the labels. Only do this if
# the image width would not end up less than 8 (this can cause problems).
imageWidth <- 8
if( ( length( assayGroupLabels ) / 2 ) > 8 ) {
	imageWidth <- length( assayGroupLabels ) / 2
}

# Stater value for overall image height,
imageHeight <- 8
# and image margin height.
marginHeight <- 8

# See how long the longest assay group label is. If it's over 16 characters,
# need to make the image height and margin height larger.
# Get the lengths of all the assay group labels.
assayGroupLabelLengths <- sapply( assayGroupLabels, function( x ) nchar( x ) )
# How long is the longest one?
longestLabel <- assayGroupLabelLengths[ which.max( assayGroupLabelLengths ) ]

# Some messing around with image and margin height to get the column (assay
# group) labels to fit on the page. This is not ideal and in the long run we
# should probably think harder about a better solution.
if( longestLabel / 3 > 8 ) {
	imageHeight <- ( longestLabel / 3 )
	marginHeight <- ( longestLabel / 3 )
}

# Make the heatmap.
message( paste( "Drawing heatmap in", heatmapFilename ) )
pdf( heatmapFilename, height=imageHeight, width=imageWidth )
heatmap.2( 
	as.matrix( top100geneFPKMs ), 
	col = colours, 
	labRow = top100geneNames,
	labCol = assayGroupLabels,
	key = FALSE,
	trace = "none",
	cexRow = 0.4,
	cexCol = 0.7,	# hardcoding for now, may need to make this dynamic but requires thinking about.
	margins = c( marginHeight, 6 )
)
invisible( dev.off() )

