#!/usr/bin/env Rscript

# Script to read in a matrix of FPKMs, run quantile normalization, 
# and write out the normalized matrix to a new file. 

# Load limma package.
suppressMessages( library( limma ) )

# Function check_file_exists -- quite if the file doesn't exist.
check_file_exists <- function( filename ) {

    if( !file.exists( filename ) ) {
        stop( paste(
            "Cannot file:",
            filename
        ) )
    }
}

# Get command line arguments -- filename for input and filename for output.
args <- commandArgs( TRUE )

# Make sure we have the correct number of arguments.
if( length( args ) != 3 ) {

    stop( "Usage\n\tgxa_quantileNormalization.R <input matrix file> <assay groups to column headings file> <output matrix file>" )

} else {

    inputFile <- args[ 1 ]
    assayGroupsMappingsFile <- args[ 2 ]
    outputFile <- args[ 3 ]
}

# Check that the input file exists.
check_file_exists( inputFile )

# Read in the matrix file and add rownames.
fpkmsDF <- read.delim( inputFile, header = TRUE )
rownames( fpkmsDF ) <- fpkmsDF[ , 1 ]
fpkmsDF[ , 1 ] <- NULL

# Read in the mappings between assay groups and column headings.
assayGroupMappingsDF <- read.delim( assayGroupsMappingsFile, header = TRUE, stringsAsFactors = FALSE )

# Find the unique assay groups.
assayGroups <- unique( assayGroupMappingsDF$AssayGroupID )
names( assayGroups ) <- assayGroups

# Create a list mapping assay groups to vectors of column headings.
assayGroupMappingsList <- lapply( assayGroups, function( assayGroup ) {

    # Get the row indices of this assay group in the data frame.
    rowIndices <- which( assayGroupMappingsDF$AssayGroupID == assayGroup )

    # Return the column headings at those row indices.
    make.names( assayGroupMappingsDF$ColumnHeading[ rowIndices ] )
})

# Create a list of data frames containing normalized data, one data frame for
# each assay group.
normalizedFpkmsList <- lapply( assayGroupMappingsList, function( columnHeadings ) {
    
    # Get the data frame columns for this assay group. Use drop = FALSE to catch
    # single-assay assay groups.
    thisAssayGroupColumns <- fpkmsDF[ , columnHeadings, drop = FALSE ]

    # Convert the dataframe to a matrix prior to normalization
    thisAssayGroupMatrix <- as.matrix( thisAssayGroupColumns )

    # Run quantile normalization.
    thisAssayGroupNormalizedMatrix <- normalizeBetweenArrays( thisAssayGroupMatrix, method = "quantile" )

    # Return a dataframe with the gene IDs in a column using the normalized matrix.
    data.frame( thisAssayGroupNormalizedMatrix )
})

# Combine the data frames into one data frame.
normalizedFpkmDF <- data.frame( normalizedFpkmsList[ names( normalizedFpkmsList ) ] )
# Remove the assay group IDs from the column headings.
colnames( normalizedFpkmDF ) <- sub( "^g\\d+\\.", "", colnames( normalizedFpkmDF ) )

# Add the gene IDs as the first column.
normalizedFpkmDF <- data.frame( GeneID = rownames( normalizedFpkmDF ), normalizedFpkmDF )

# Write the data frame containing normalized data to a new file.
write.table( normalizedFpkmDF, file = outputFile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE )