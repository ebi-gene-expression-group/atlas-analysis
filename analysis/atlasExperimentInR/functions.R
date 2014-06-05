source("/ebi/microarray/home/mkeays/Atlas/git/atlasprod/analysis/atlasExperimentInR/Analytics.R")

# Required packages.
library( XML )
library( plyr )
library( Biobase )
library( GenomicRanges )

# Some bits of config.
# ArrayExpress load directories -- where SDRFs live.
ae2experiments <- "/ebi/microarray/home/arrayexpress/ae2_production/data/EXPERIMENT";

# Directory with gene annotations.
bioentityPropertiesEnsemblDir <- "/ebi/microarray/home/atlas3-production/bioentity_properties/ensembl";

# summarizeAtlasExperiment
# 	- Main function for the package. Takes an experiment accession and a
# 	directory path where it can find Atlas XML config and expressions matrices.
# 	- Returns a list of ExpressionSet and/or MAList and/or SummarizedExperiment objects.
summarizeAtlasExperiment <- function( experimentAccession, atlasExperimentDirectory ) {
	
	# Atlas XML config file name.
	atlasExperimentXMLfile <- paste( experimentAccession, "-configuration.xml", sep="" )
	atlasExperimentXMLfile <- file.path( atlasExperimentDirectory, atlasExperimentXMLfile )

	# Die if we can't find the XML config file.
	if( !file.exists( atlasExperimentXMLfile ) ) {
		stop( paste( "XML config file", atlasExperimentXMLfile, "does not exist or is not readable." ) )
	}

	# Parse the XML file.
	experimentXMLlist <- parseAtlasXML( atlasExperimentXMLfile )

	# Get the pipeline code from the experiment accession e.g. MTAB, MEXP
	pipeline <- gsub( "E-", "", experimentAccession )
	pipeline <- gsub( "-\\d+", "", pipeline )

	# Filename for SDRF.
	sdrfBasename <- paste( experimentAccession, ".sdrf.txt", sep="" )
	
	# Complete path to SDRF file.
	sdrfPath <- file.path( ae2experiments, pipeline, experimentAccession, sdrfBasename )

	# Check SDRF exists, die if not.
	if( !file.exists( sdrfPath ) ) {
		stop( paste( "SDRF", sdrfPath, "does not exist or is not readable." ) )
	}

	# Get the experiment type from the parsed XML.
	atlasExperimentType <- experimentXMLlist$experimentType

	# Parse the SDRF.
	atlasSDRF <- parseSDRF( sdrfPath, atlasExperimentType )

	# Get the list of Analytics objects.
	allAnalytics <- experimentXMLlist$allAnalytics

	# Create path to analysis methods file.
	analysisMethodsFile <- paste( experimentAccession, "-analysis-methods.tsv", sep="" )
	analysisMethodsFile <- file.path( atlasExperimentDirectory, analysisMethodsFile )

	# Next step is to go through the analytics objects created from the XML,
	# pull out the right rows form the SDRF, get the right expressions matrix,
	# the gene annotations, and make the Bioconductor object (ExpressionSet,
	# MAList, or SummarizedExperiment).
	atlasExperimentSummary <- lapply( allAnalytics, function( analytics ) {
	
		analyticsSDRF <- createAnalyticsSDRF( analytics, atlasSDRF )

		# Get the expressions.
		expressionsDF <- getExpressions( analytics, atlasExperimentType, experimentAccession, atlasExperimentDirectory )
		
		# If this is a one-colour microarray experiment, make an ExpressionSet.
		if( grepl( "microarray_1colour", atlasExperimentType ) ) {
			
			# Create an ExpressionSet.
			expressionSet <- createExpressionSet( expressionsDF, analyticsSDRF, analysisMethodsFile )
			
			# Return it.
			return( expressionSet )
		}

		# If this is an RNA-seq experiment, make a SummarizedExperiment.
		else if( grepl( "rnaseq", atlasExperimentType ) ) {
			
			# Create a SummarizedExperiment.
			summarizedExperiment <- createSummarizedExperiment( expressionsDF, analyticsSDRF, analysisMethodsFile )
		}
			
	})

	return( atlasExperimentSummary )

}


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

	# If this is a baseline RNA-seq experiment there won't be any raw counts,
	# so can't create a BioC object.
	if( grepl( "^rnaseq\\w*baseline$", atlasExperimentType ) ) {
		stop( paste( "Cannot create R objects for experiment type \"", atlasExperimentType, "\".", sep="" ) )
	}
	
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

			assayNameColIndex <- grep( "Hybridi[sz]ation\\s?Name", completeSDRF[ 1, ] )

			# Check that we got something.
			if( length( assayNameColIndex ) != 1 ) {
				stop( "Did not find an Assay Name column or a Hybridization Name column, cannot continue." )
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
	if( length( duplicateColIndices ) > 0 ) {
		subsetSDRF <- subsetSDRF[ , -duplicateColIndices ]
		newColNames <- newColNames[ -duplicateColIndices ]
	}
	
	# Remove the first row of the SDRF (this is the old column headings)
	subsetSDRF <- subsetSDRF[ -1, ]

	# Add the new column names as the column headings.
	colnames( subsetSDRF ) <- newColNames
	
	# Remove duplicated rows, which occur e.g. if an assay has more than one file.
	duplicateRowIndices <- which( duplicated( subsetSDRF ) )
	if( length( duplicateRowIndices ) > 0 ) {
		subsetSDRF <- subsetSDRF[ -duplicateRowIndices, ]
	}

	# Make assay names "R-safe".
	subsetSDRF$AssayName <- make.names( subsetSDRF$AssayName )

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


# getExpressions
# 	- Takes an Analytics object, experiment type, experiment accession, and
# 	path to directory containing expressions file.
# 	- Returns a data frame containing the expressions.
getExpressions <- function( analytics, atlasExperimentType, experimentAccession, atlasExperimentDirectory ) {
	
	
	# Is this a 1-colour microarray experiment?
	if( grepl( "microarray_1colour", atlasExperimentType ) ) {
		
		# Need the array design as it makes up the file name.
		arrayDesign <- platform( analytics )
		
		# Create the file name of the normalized expressions.
		expressionsFile <- paste( experimentAccession, "_", arrayDesign, "-normalized-expressions.tsv.undecorated", sep="" )
		
		# Add the full path to the file.
		expressionsFile <- file.path( atlasExperimentDirectory, expressionsFile )

		# Check that the expressions file exists, die if not.
		if( !file.exists( expressionsFile ) ) {
			stop( paste( "Expressions file \"", expressionsFile, "\" does not exist. Cannot continue.", sep="" ) )
		}
		
		# Read the expressions file.
		expressionsDF <- read.delim( expressionsFile, header=TRUE, stringsAsFactors=FALSE )
		
		# Add design elements as row names.
		rownames( expressionsDF ) <- expressionsDF[,1]
		
		# Remove design elements column as not needed.
		expressionsDF[,1] <- NULL
	}
	# Is this an RNA-seq experiment?
	else if( grepl( "rnaseq", atlasExperimentType ) ) {
		
		# Create the file name of the raw counts.
		expressionsFile <- paste( experimentAccession, "-raw-counts.tsv.undecorated", sep="" )
		
		# Add the full path to the file.
		expressionsFile <- file.path( atlasExperimentDirectory, expressionsFile )
		
		# Check that the expressions file exists, die if not.
		if( !file.exists( expressionsFile ) ) {
			stop( paste( "Expressions file \"", expressionsFile, "\" does not exist. Cannot continue.", sep="" ) )
		}
		
		# Read the expressions file.
		expressionsDF <- read.delim( expressionsFile, header=TRUE, stringsAsFactors=FALSE )
		
		# Add gene IDs as row names.
		rownames( expressionsDF ) <- expressionsDF[,1]
		
		# Remove now unwanted column.
		expressionsDF[,1] <- NULL
	}
	# Otherwise, don't recognise this experiment type so die.
	else {
		stop( paste( "Don't know how handle experiment type \"", atlasExperimentType, "\". Cannot continue.", sep="" ) )
	}
	
	# Return the data frame with expressions.
	return( expressionsDF )
}


# createAnalyticsSDRF
# 	- Take an Analytics object and a parsed SDRF.
# 	- Return a new SDRF data frame with just the assays from the Analytics
# 	object, as well as a column denoting which assay group each assay belongs
# 	to.
createAnalyticsSDRF <- function( analytics, atlasSDRF ) {

	# Get the assay groups.
	assayGroups <- assay_groups( analytics )

	# Get the SDRF rows for these assay groups.
	assayGroupSDRFs <- lapply( assayGroups, function( assayGroup ) {
		
		# Get the assay names.
		assayNames <- assays( assayGroup )

		# TODO: strip .Cy* from 2-colour assay names? Check diffAtlas_DE_limma.R

		# Get the SDRF rows for these assays.
		atlasSDRF[ which( atlasSDRF$AssayName %in% assayNames ), ]
	})

	# Make a new data frame from the SDRF chunks list.
	analyticsSDRF <- ldply( assayGroupSDRFs, data.frame )

	# Change the name of the first column.
	colnames( analyticsSDRF )[1] <- "AtlasAssayGroup"

	# Make the assay names the row names.
	rownames( analyticsSDRF ) <- analyticsSDRF$AssayName
	analyticsSDRF$AssayName <- NULL

	# Sort the rows by assay name.
	analyticsSDRF <- analyticsSDRF[ sort( rownames( analyticsSDRF ) ) , ]

	return( analyticsSDRF )
}


# createExpressionSet
# 	- Take a data frame of normalized expressions and a parsed SDRF.
# 	- Return an ExpressionSet object.
createExpressionSet <- function( expressionsDF, analyticsSDRF, analysisMethodsFile ) {

	# Only select columns with assay names in our SDRF -- these are the
	# ones that passed QC and are still in the XML.
	expressionsDF <- expressionsDF[ , rownames( analyticsSDRF ) ]

	# Turn data frame into matrix.
	expressionsMatrix <- as.matrix( expressionsDF )
	
	# Create a new AssayData object
	expressionData <- assayDataNew( storage.mode = "lockedEnvironment", exprs = expressionsMatrix )
	
	# Add feature names (probe set names).
	featureNames( expressionData ) <- rownames( expressionsMatrix )
	
	# Add sample names (assay names).
	sampleNames( expressionData ) <- colnames( expressionsMatrix )
	
	# Add the SDRF data.
	phenoData <- new( "AnnotatedDataFrame", data = analyticsSDRF )
	
	# Add featureData -- this is just the probe set names for now.
	featureData <- new( "AnnotatedDataFrame", data = data.frame( probeSets = rownames( expressionsMatrix ) ) )
	featureNames( featureData ) <- rownames( expressionsMatrix )
	
	# Add analysis methods.
	analysisMethodsList <- readArrayAnalysisMethods( analysisMethodsFile )
	# Add this to a MIAME object.
	exptData <- new( "MIAME", preprocessing = analysisMethodsList )

	# Create new ExpressionSet.
	return( new( "ExpressionSet", 
		assayData = expressionData, 
		phenoData = phenoData, 
		featureData = featureData, 
		experimentData = exptData 
	) )
}


createSummarizedExperiment <- function( expressionsDF, analyticsSDRF, analysisMethodsFile ) {
	
	# Only select columns with assay names in our SDRF.
	expressionsDF <- expressionsDF[ , rownames( analyticsSDRF ) ]

	# Turn data frame into matrix.
	expressionsMatrix <- as.matrix( expressionsDF )

	# Turn SDRF into a DataFrame.
	analyticsSDRF <- DataFrame( analyticsSDRF )

	# Get analysis methods
	analysisMethodsList <- readSeqAnalysisMethods( analysisMethodsFile )

	# Create SummarizedExperiment
	summarizedExperiment <- SummarizedExperiment( assays = expressionsMatrix, colData = analyticsSDRF, exptData = analysisMethodsList )

	# Return it.
	return( summarizedExperiment )
}


readArrayAnalysisMethods <- function( analysisMethodsFile ) {
	
	# Read the analysis methods file.
	analysisMethodsDF <- read.delim( analysisMethodsFile, header=FALSE, stringsAsFactors=FALSE )
	
	# Find the row that has the normalization method.
	normalizationRow <- which( analysisMethodsDF[ , 1 ] == "Normalization" )
	
	# Get the text on this row in the second column.
	normalizationText <- analysisMethodsDF[ normalizationRow, 2 ]
	
	# Reorganise text.
	normalizationText <- sub( "^(.*) <a href=(.*)>(.*)</a> (version.*). <a.*", "\\1 \\3 (\\2) \\4", normalizationText)
	
	# Make a list containing the normalization text.
	analysisMethodsList <- list( normalization = normalizationText )
	
	# Return it.
	return( analysisMethodsList )
}


readSeqAnalysisMethods <- function( analysisMethodsFile ) {
	
	# Read the analysis methods file.
	analysisMethodsDF <- read.delim( analysisMethodsFile, header=FALSE, stringsAsFactors=FALSE )
	
	# Get the row with iRAP information.
	irapRow <- grep( "Pipeline version", analysisMethodsDF[ , 1 ] )

	# Get the text for iRAP.
	irapInfo <- analysisMethodsDF[ irapRow, 2 ]

	# Reorganise iRAP text.
	irapInfo <- sub( "^<a href=(.*)>(.*)</a> (.*) ", "\\2 version \\3 (\\1)", irapInfo, )

	# Get the rows with filtering info.
	filteringRows <- grep( "Filtering Step", analysisMethodsDF[ , 1 ] )

	# Get the text for the filtering steps.
	filteringInfo <- analysisMethodsDF[ filteringRows, 2 ]

	# Get the row woth mapping info.
	mappingRow <- grep( "Read Mapping", analysisMethodsDF[ , 1 ] )

	# Get the text for the mapping info.
	mappingInfo <- analysisMethodsDF[ mappingRow, 2 ]

	# Get the row with quantification info.
	quantRow <- grep( "Quantification", analysisMethodsDF[ , 1 ] )
	
	# Get the quantification info.
	quantInfo <- analysisMethodsDF[ quantRow, 2 ]

	analysisMethodsList <- SimpleList( 
		pipeline = irapInfo,
		filtering = filteringInfo,
		mapping = mappingInfo,
		quantification = quantInfo
	)
	
	return( analysisMethodsList )
}

















