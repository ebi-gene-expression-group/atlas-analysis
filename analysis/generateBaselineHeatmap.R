#!/usr/bin/env Rscript
#
# Script to create a heatmap for a baseline Expression Atlas experiment.

# Load ExpressionAtlasInternal package.
suppressMessages( library( ExpressionAtlasInternal ) )
# Load genefilter package to get rowVars function.
suppressMessages( library( genefilter ) )
# Load gplots package for heatmap.2 function.
suppressMessages( library( gplots ) )
# Load RColorBrewer package for nice colours.
suppressMessages( library( RColorBrewer ) )

suppressMessages( library(optparse))

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

get_median_expressions <- function( dataFrame ) {
    startCol <- 3
    expressionLevels <- as.matrix(dataFrame[ , startCol:ncol( dataFrame ) ])
    if( !all( apply( expressionLevels, c(1,2), function( x ) { grepl( ",", x ) } ) ) ) {
        stop( "your values are not comma-separated lists of quartiles. Please check." )
    }
    medians <- apply(expressionLevels, c(1,2), function(x){suppressWarnings(as.numeric(strsplit( x, "," )[[1]][3]))})
    return (
		data.frame( dataFrame[ , 1:(startCol - 1), drop=FALSE ], medians, stringsAsFactors=FALSE )
	)
}

###############################
# Script start.

args <- parse_args(OptionParser(option_list= list(
	make_option(
		c("-i", "--input"),
		help="Input tsv file : Ensembl identifier, gene name, data columns with quartiles"
	), #fpkmsMatrixFile
	make_option(
	 c("-c", "--configuration"),
	 help="Configuration file"
	 ), #experimentConfigFile
	 make_option(
		 c("-o", "--output"),
		 help="Where to save the output PDF file"
	 )
)))
check_file_exists(args$input)
check_file_exists(args$configuration)


# Parse XML config to a list.
experimentConfigList <- parseAtlasConfig( args$configuration )

# Get the AssayGroup objects from the experiment config ready to use later.
# First get list of analytics.
experimentAnalytics <- experimentConfigList$allAnalytics
# Now get the RNA-seq analytics list (there should be the only analytics element).
rnaseqAnalytics <- experimentAnalytics$rnaseq
# Now get the list of RNA-seq assay groups.
rnaseqAssayGroups <- assay_groups( rnaseqAnalytics )

# Check that all assay groups have a label, we cannot continue without one.
invisible( lapply( rnaseqAssayGroups, function( assayGroup) {
	if( length( assay_group_label( assayGroup ) ) == 0 ) {
		stop( paste( "Assay group", assay_group_id( assayGroup ), "does not have a label. Cannot continue." ) )
	}
} ) )
#Expected format:
#gene id | gene name | g1 | .. gn
dataFrame<-read.delim( args$input, stringsAsFactors = FALSE, header = TRUE )

expressionsNumeric <- get_median_expressions(dataFrame)
# Assign gene IDs as row names.
rownames( expressionsNumeric ) <- expressionsNumeric[[1]]
# Remove the gene id column.
expressionsNumeric[[1]] <- NULL

# Create data frame of just Ensembl IDs and gene names. We have to use Ensembl
# IDs as row names in R because it doesn't allow non-unique row names.
geneIDsToGeneNames <- data.frame( id = expressionsNumeric[[1]], stringsAsFactors=FALSE )
# Add the Ensembl gene IDs as the row names.
rownames( geneIDsToGeneNames ) <- rownames( expressionsNumeric )

# Replace any empty gene names with the gene ID.
emptyGeneNameIndices <- which( geneIDsToGeneNames == "" )
geneIDsToGeneNames[ emptyGeneNameIndices , ] <- rownames( geneIDsToGeneNames )[ emptyGeneNameIndices ]

# Now remove the Gene.Name column from the data frame.
expressionsNumeric[[1]] <- NULL

# To create the heatmap, we need to take the top 100 most variable genes. To do
# this, we will use the "rowVars" function from the genefilter package, which
# calculates the variance of each row of a data frame.
if(dim(expressionsNumeric)[2] > 1) {
	rowVariances <- rowVars( expressionsNumeric )
} else {
	# can't calculate variance of one row, use the FPKM values (alhough the heatmap isn't as valuable then)
	rowVariances <- expressionsNumeric[1]
}

# Get the indices of the top 100 most variable genes.
chosenColumnIndices<- order( rowVariances, decreasing = TRUE )[1:100]
top100geneExpressions <- expressionsNumeric[ chosenColumnIndices , ]

# Get the gene names for the top 100 gene IDs, to use as labels for the
# heatmape rows.
top100geneNames <- geneIDsToGeneNames[ chosenColumnIndices, ]

# Scale and center the expression levels using Z-score transformation, so they
# have mean 0 and standard deviation 1. This makes the clustering and heatmap
# colours look better than leaving them as they are (very long tail).
top100geneExpressions <- t( scale( t( top100geneExpressions )))

# Get the assay group labels to use as the labels for the heatmap columns.
assayGroupLabels <- lapply( colnames( top100geneExpressions ), function( assayGroupID ) {
	assayGroup <- rnaseqAssayGroups[[ assayGroupID ]]
	assayGroupLabel <- assay_group_label( assayGroup )
} )

# Some nice colours.
colours <- colorRampPalette( brewer.pal( 9, "Blues" ) )( 100 )

# Make the heatmap width a function of the number of assay groups. If the
# columns are too narrow it's impossible to read the labels. Only do this if
# the image width would not end up less than 8 (this can cause problems).
imageWidth <- 8
if( ( length( assayGroupLabels ) / 2 ) > 8 ) {
	imageWidth <- length( assayGroupLabels ) / 2
}

# See how long the longest assay group label is. If it's over 16 characters,
# need to make the image height and margin height larger.
# Get the lengths of all the assay group labels.
# How long is the longest one?
longestLabel <- max(unlist(lapply( assayGroupLabels, function( x ) nchar( x ) )))

# Some messing around with image and margin height to get the column (assay
# group) labels to fit on the page. This is not ideal and in the long run we
# should probably think harder about a better solution.
if( longestLabel / 3 > 8 ) {
	imageHeight <- ( longestLabel / 3 )
	marginHeight <- ( longestLabel / 3 )
} else {
	imageHeight <- 8
	marginHeight <- 8
}

# Make the heatmap.
pdf( args$output, height=imageHeight, width=imageWidth )
heatmap.2(
	as.matrix( top100geneExpressions ),
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
