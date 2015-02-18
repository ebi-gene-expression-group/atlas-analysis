#!/usr/bin/env Rscript
#
# Script to create a heatmap for a baseline Expression Atlas experiment.

# Load atlasSiteConfig package to read YAML site config.
suppressMessages( library( atlasSiteConfig ) )
# Load atlasConfigParser package to read XML config.
suppressMessages( library( atlasConfigParser ) )
# Load genefilter package to get rowVars function.
suppressMessages( library( genefilter ) )
# Load gplots package for heatmap.2 function.
suppressMessages( library( gplots ) )
# Load RColorBrewer package for nice colours.
suppressMessages( library( RColorBrewer ) )

atlasProdDir <- Sys.getenv( "ATLAS_PROD" )

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

get_ensgene_filename <- function( species, atlasProdDir ) {
	
	# Parse site config to a list.
	siteConfig <- createAtlasSiteConfig()

	# Get the Ensembl bioentity properties directory.
	ensemblDirectory <- siteConfig$bioentity_properties_ensembl

	# Create path to ensgene file.
	ensgeneFilePath <- file.path( atlasProdDir, ensemblDirectory, paste( species, "ensgene", "tsv", sep="." ) )

	return( ensgeneFilePath )
}

make_species_specific_data_frame <- function( speciesFPKMsFile, speciesEnsgeneFile ) {

	# Read in the files now we have them.
	message( paste( "Reading FPKMs from", speciesFPKMsFile, "..." ) )
	fpkmsDataFrame <- read.delim( speciesFPKMsFile, stringsAsFactors=FALSE, header=TRUE )
	message( "Successfully read FPKMs" )

	message( paste( "Reading Ensembl gene annotations from", speciesEnsgeneFile, "..." ) )
	ensgene <- read.delim( speciesEnsgeneFile, stringsAsFactors=FALSE, header=TRUE )
	message( "Successfully read Ensembl gene annotations" )

	message( "Adding gene names to FPKMs data frame..." )

	# Just need the "ensgene" and "symbol" columns from the bioentity
	# properties file. These are the gene names and gene IDs, respectively.
	ensgene <- ensgene[ , c( "ensgene", "symbol" ) ]

	# Get all the gene names from ensgene in the same order as the gene IDs in
	# the FPKMs matrix. If there is no name for a gene, use an emtpy string.
	geneNames <- sapply( fpkmsDataFrame$GeneID, function( geneID ) {
		
		# Get the index of the gene ID in the ensgene data frame.
		ensGeneIdx <- which( ensgene$ensgene == geneID )
		
		# If the index has length >0 then there was a match.
		if( length( ensGeneIdx ) > 0 ) {
			# Get the gene name at that index.
			name <- ensgene$symbol[ ensGeneIdx ]
		} else {
			# If there wasn't a match, use an empty string.
			name <- ""
		}
	})

	# Add the gene names to the FPKMs data frame.
	fpkmsDataFrame$Gene.Name <- geneNames
	
	# Need to do some reshuffling and renaming of columns in the FPKMs data frame.
	# Count the columns.
	numCol <- ncol( fpkmsDataFrame )
	# Get column indices for GeneID and Gene.Name columns
	geneIDcol <- which( colnames( fpkmsDataFrame ) == "GeneID" )
	geneNameCol <- which( colnames( fpkmsDataFrame ) == "Gene.Name" )

	# The assay group columns are the rest of the columns.
	assayGroupCols <- 1:numCol
	assayGroupCols <- assayGroupCols[ -c( geneIDcol, geneNameCol ) ]
	
	# Re-shuffle the column order.
	fpkmsDataFrame <- fpkmsDataFrame[ , c( geneIDcol, geneNameCol, assayGroupCols ) ]

	# Rename the GeneID column so it's the same as the decorated FPKMs matrices
	# for single-species experiments when it's read in.
	colnames( fpkmsDataFrame )[ 1 ] <- "Gene.ID"

	message( "Successfully added gene names to FPKMs data frame" )

	return( fpkmsDataFrame )
}

###############################
# Script start.


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
check_file_exists( atlasExperimentDirectory )

# Build the XML config file name and the FPKM matrix filename.
# Need to get the experiment accession out of the path passed.
experimentAccession <- basename( atlasExperimentDirectory )

# Log accession.
message( paste( "Experiment accession is:", experimentAccession ) )

# XML config file.
experimentConfigFile <- file.path( atlasExperimentDirectory, paste( experimentAccession, "-configuration.xml", sep="" ) )
# Check it exists.
check_file_exists( experimentConfigFile )

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



# If species was provided, we want to append to accession, but there won't
# always be a species provided, so we have a new variable that can be either
# just accession, or accession_species. This is used for FPKM matrix and
# heatmap filenames.
experimentAccessionForFilename <- experimentAccession

# If species was provided...
if( exists( "species" ) ) {
	# Log species.
	message( paste( "Species is:", species ) )
	
	# Add the species to the accession, for filenames.
	experimentAccessionForFilename <- paste( experimentAccessionForFilename, species, sep="_" )
	
	# The undecorated FPKM matrix for this species.
	speciesFPKMsFile <- file.path( atlasExperimentDirectory, paste( experimentAccessionForFilename, ".tsv.undecorated", sep="" ) )

	# Check it exists.
	check_file_exists( speciesFPKMsFile )
	
	# Get the Ensembl bioentity properties "ensgene" file name.
	speciesEnsgeneFile <- get_ensgene_filename( species, atlasProdDir )
	check_file_exists( speciesEnsgeneFile )

	fpkmsDataFrame <- make_species_specific_data_frame( speciesFPKMsFile, speciesEnsgeneFile )

} else {
	# FPKMs matrix file.
	fpkmsMatrixFile <- file.path( atlasExperimentDirectory, paste( experimentAccession, ".tsv", sep="" ) )

	# Check the FPKMs matrix exists.
	check_file_exists( fpkmsMatrixFile )

	# Read in the FPKMs.
	message( paste( "Reading FPKMs from", fpkmsMatrixFile, "..." ) )
	fpkmsDataFrame <- read.delim( fpkmsMatrixFile, stringsAsFactors=FALSE, header=TRUE )
	message( "Successfully read FPKMs" )
}

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

# Check that there is more than one column of FPKMs. We can't make a heatmap
# with less than two columns of data.
if( ncol( fpkmsDataFrame ) < 2 ) {
	# Warn what has happened.
	warning( "Less than two columns of FPKMs found. Cannot continue." )
	
	# Exit without exit code.
	q( save="no" )
}

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
heatmapFilename <- paste( experimentAccessionForFilename, "-heatmap.pdf", sep="" )
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



