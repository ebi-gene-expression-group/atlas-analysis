source("/ebi/microarray/home/mkeays/Atlas/git/atlasprod/analysis/atlasExperimentInR/Analytics.R")

# parseAtlasXML
# 	- Read Atlas XML config file and return a list of Analytics objects, as
# 	well as the experiment type.
# 	- The function returns a list with two elements: the list of Analytics
# 	objects, and the experiment type from the XML.
parseAtlasXML <- function(filename) {
	
	# Load XML package.
	library( XML )

	# Read the XML file.
	xmlTree <- xmlInternalTreeParse( filename )

	# Get the configuration node -- the root of the tree.
	configNode <- xmlRoot( xmlTree )

	# Get the config node attributes.
	configNodeAttrs <- xmlAttrs( configNode )

	# Get the Atlas experiment type. We'll use this later to decide whether to
	# create an ExpressionSet/MAList or SummarizedExperiment.
	atlasExperimentType <- configNodeAttrs[[ "experimentType" ]]
	
	# Go through the analytics node(s).
	# Get them from the configuration node.
	allAnalyticsNodes <- xmlElementsByTagName( configNode, "analytics" )

	# Create a list of Analytics objects.
	allAnalytics <- lapply( allAnalyticsNodes, function( analyticsNode ) {
		analyticsObject <- new( "Analytics", atlasExperimentType, analyticsNode )
	})
	
	names( allAnalytics ) <- sapply( allAnalytics,
		function( analytics ) {
			platform( analytics )
		}
	)

	# Since we can't return more than one thing, make a list to put the
	# analytics list in as well as the experiment type.
	parsedXMLlist <- list()
	
	parsedXMLlist$allAnalytics <- allAnalytics
	parsedXMLlist$experimentType <- atlasExperimentType

	return( parsedXMLlist )
}


# parseSDRF
# 	- Take an SDRF filename and the experiment type from the XML config, and
# 	return a subset of the SDRF containing only the assay names (or ENA runs),
# 	the Characteristics, and FactorValue columns. 
# 	- It tries to include unit columns as well.
# 	- It renames the columns to remove e.g. Characteristics[].
# 	- It removes duplicated columns, e.g. if genotype is a Characteristic and a
# 	Factor.
# 	- It returns the new "SDRF" as a data frame.
parseSDRF <- function( filename, atlasExperimentType ) {

	# Read in the SDRF file. Set header=FALSE because we don't want the column
	# headings to be made "R-safe" -- this confuses things when we're trying to
	# find the Charactersitic and Factor names. Set stringsAsFactors=FALSE so
	# that we can use grep.
	completeSDRF <- read.delim( filename, header = FALSE, stringsAsFactors = FALSE )

	# Get the Characteristics column indices, and any unit columns next to them.
	charColIndices <- grep( "Characteristics", completeSDRF[ 1, ] )
	charColIndices <- addUnitCols( charColIndices, completeSDRF )
	
	# Get the Factor column indices, an any unit columns next to them.
	factorColIndices <- grep( "Factor\\s?Value", completeSDRF[ 1, ] )
	factorColIndices <- addUnitCols( factorColIndices, completeSDRF )
	
	# Get the column index for assay names. For microarray data, this is "Assay
	# Name" or "Hybridization Name". For RNA-seq data, this is "Comment[ENA_RUN]"
	if( grepl( "rnaseq", atlasExperimentType ) ) {
		
		assayNameColIndex <- grep( "Comment\\s?\\[\\s?ENA_RUN\\s?\\]", completeSDRF[ 1, ] )

		if( length( assayNameColIndex ) != 1 ) {
			stop( "Did not find Comment[ ENA_RUN ] column in SDRF." )
		}
	}
	else {
		
		assayNameColIndex <- grep( "Assay\\s?Name", completeSDRF[ 1, ] )

		if( length( assayNameColIndex ) != 1 ) {
			print( "Did not find Assay Name column in SDRF. Looking for Hybridization Name instead ..." )

			assayNameColIndex <- grep( "Hybridi[sz]ation\\s?Name", completeSDRF[ 1, ] )

			# Check that we got something.
			if( length( assayNameColIndex ) != 1 ) {
				stop( "Did not find Hybridization Name column either, cannot continue." )
			}
		}

		# For two-colour array data, also want to get the label column.
		if( grepl( "2colour", atlasExperimentType ) ) {

			labelColIndex <- which( completeSDRF[ 1, ] == "Label" )
		}
	}
	
	# Now we should have everything we need to get the right columns and make a
	# more friendly SDRF.
	if( grepl( "2colour", atlasExperimentType ) ) {

		subsetSDRF <- completeSDRF[ , c( assayNameColIndex, labelColIndex, charColIndices, factorColIndices ) ]
	}
	else {
		
		subsetSDRF <- completeSDRF[ , c( assayNameColIndex, charColIndices, factorColIndices ) ]
	}
	
	# Next thing is to name the columns so they have nice names.
	newColNames <- gsub( "Characteristics\\s?\\[", "", subsetSDRF[1,] )
	newColNames <- gsub( "Factor\\s?Value\\s?\\[", "", newColNames )
	newColNames <- gsub( "Unit\\s?\\[", "", newColNames )
	newColNames <- gsub( "\\s?\\]", "", newColNames )
	newColNames[ 1 ] <- "AssayName"
	
	# Replace spaces with underscores.
	newColNames <- gsub( " ", "_", newColNames )

	# Now we've got the new names for the columns, check if any are the same
	# (use "tolower" function to convert all to lower case).
	duplicateColIndices <- which( duplicated( tolower( newColNames ) ) )

	# Remove the duplicated columns from the new column names and the subset
	# SDRF.
	subsetSDRF <- subsetSDRF[ , -duplicateColIndices ]
	newColNames <- newColNames[ -duplicateColIndices ]
	
	# Remove the first row of the SDRF (this is the old column headings)
	subsetSDRF <- subsetSDRF[ -1, ]

	# Add the new column names as the column headings.
	colnames( subsetSDRF ) <- newColNames

	# For 2-colour experiments, merge the AssayName and Label columns.
	if( grepl( "2colour", atlasExperimentType ) ) {

		subsetSDRF$AssayName <- paste( subsetSDRF$AssayName, subsetSDRF$Label, sep="." )
		
		# Remove the Label column -- the second one.
		subsetSDRF <- subsetSDRF[ , -2 ]
	}
	
	# Remove duplicated rows, which occur e.g. if an assay has more than one file.
	duplicateRowIndices <- which( duplicated( subsetSDRF ) )
	subsetSDRF <- subsetSDRF[ -duplicateRowIndices, ]

	# Return the subset SDRF.
	return( subsetSDRF )
}


# addUnitCols
# 	- Given a vector of column indices and a data frame with the complete SDRF,
# 	return a vector of column indices containing the original ones plus any
# 	Unit[] columns that are next to them.
addUnitCols <- function( colIndices, SDRF ) {

	# Get the indices of unit columns. 
	unitCols <- unlist(
		
		# Go through each column index provided.
		sapply(
			colIndices,
			function( colNumber ) {
				
				# Look at the column to the right of it.
				nextCol = colNumber + 1

				# Only try if this is not the very last column.
				if( nextCol <= ncol( SDRF ) ) {

					# If it's a unit column, return the index.
					if( grepl( "Unit", SDRF[ 1, nextCol ] ) ) {
						nextCol
					}
				}
			}
		)
	)
	
	# Combine the vector of unit column indices with the original one, and then
	# sort them so the unit column indices are next to the correct non-unit
	# column indices. E.g. if "Characteristics[ age ]" was the 2nd column out
	# of 5, that means its "Unit[ time unit ]" column is the 3rd column.
	# Combining the unit and non-unit column indices gives:
	# 	1, 2, 4, 5, 3
	# So then we sort it so it becomes:
	# 	1, 2, 3, 4, 5
	# Now when we use this vector of indices to select the Characteristics
	# columns we'll get them in the right order.
	allColIndices <- sort( c( colIndices, unitCols ) )

	return( allColIndices )
}


summarizeAtlasExperiment <- function( sdrfDataFrame, allAnalytics ) {

}
