
# parseAtlasXML
# 	- Read Atlas XML config file and return a list of Analytics objects, as
# 	well as the experiment type.
# 	- The function returns a list with two elements: the list of Analytics
# 	objects, and the experiment type from the XML.
parseAtlasXML <- function( atlasXMLconfigFile ) {
	
	# Read the XML file.
	xmlTree <- xmlInternalTreeParse( atlasXMLconfigFile )

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


deduceSiteConfigPath <- function( dbName ) {

    # We need the path of the YAML config file.
    # We can deduce this using the output of the R.home() function, which is
    # the directory where R is installed.
    rHome <- R.home()

    # Split the path on "/" characters so that we can remove everything after
    # the atlasinstall directory more easily. Use the platform's file separator
    # (in UNIX a "/").
    splitPath <- strsplit( rHome, .Platform$file.sep )[[ 1 ]]

    # Find the index of the atlasinstall_<instance> directory.
    atlasInstallIndex <- grep( "atlasinstall", splitPath )

    # Check that we got an index. If not, we must not be using an atlasinstall R.
    if( length( atlasInstallIndex ) == 0 ) {
        stop( "Could not find atlasinstall directory in RHOME path. Unable to deduce path to site config file. Please ensure you are using an atlasinstall R installation." )
    }

    # Remove everything after the atlasinstall directory from the path.
    atlasInstallSplitPath <- head( splitPath, n = atlasInstallIndex )

    # Stick the path back together again.
    atlasInstallPath <- paste( atlasInstallSplitPath, collapse=.Platform$file.sep )

	siteConfigBasename <- paste( dbName, "SiteConfig", ".yml", sep="" )
	
    # Add the path to the site config file.
    siteConfigFile <- file.path( atlasInstallPath, "atlasprod", "supporting_files", siteConfigBasename )

	return( siteConfigFile )
}


# createAtlasSiteConfig
#   - Return a list containing the Atlas site config.
createAtlasSiteConfig <- function( ) {
	
	# Get the path to the Atlas site config file.
	atlasSiteConfigFile <- deduceSiteConfigPath( "Atlas" )

    # Read in the site config file.
    atlasSiteConfig <- yaml.load_file( atlasSiteConfigFile )

    return( atlasSiteConfig )
}


# createArrayExpressSiteConfig
#	- Return a list containing the ArrayExpress site config.
createArrayExpressSiteConfig <- function( ) {
	
	# Get the path to the ArrayExpress site config file.
	arrayExpressSiteConfigFile <- deduceSiteConfigPath( "ArrayExpress" )

	# Read in the site config file
	arrayExpressSiteConfig <- yaml.load_file( arrayExpressSiteConfigFile )

	return( arrayExpressSiteConfig )
}


# summarizeAtlasExperiment
# 	- Main function for the package. Takes an experiment accession and a
# 	directory path where it can find Atlas XML config and expressions matrices.
# 	- Returns a list of ExpressionSet and/or MAList and/or SummarizedExperiment objects.
summarizeAtlasExperiment <- function( experimentAccession, atlasExperimentDirectory ) {
	
	# Quit if the accession does not look like an ArrayExpress one.
	if( !grepl( "^E-\\w{4}-\\d+$", experimentAccession ) ) {
		stop( paste( "Accession \"", experimentAccession, "\" does not look like an ArrayExpress accession. Cannot continue", sep="" ) )
	}

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
	
	# Get ArrayExpress site config to get path to experiment load directories.
	aeSiteConfig <- createArrayExpressSiteConfig( )

	ae2loadDir <- aeSiteConfig$AE2_LOAD_DIR

	# Complete path to SDRF file.
	sdrfPath <- file.path( ae2loadDir, "EXPERIMENT", pipeline, experimentAccession, sdrfBasename )

	# Check SDRF exists, die if not.
	if( !file.exists( sdrfPath ) ) {
		stop( paste( "SDRF", sdrfPath, "does not exist or is not readable." ) )
	}

	# Get the experiment type from the parsed XML.
	atlasExperimentType <- experimentXMLlist$experimentType

	# Parse the SDRF.
	atlasSDRF <- .parseSDRF( sdrfPath, atlasExperimentType )

	# Get the list of Analytics objects.
	allAnalytics <- experimentXMLlist$allAnalytics

	# Next step is to go through the analytics objects created from the XML,
	# pull out the right rows form the SDRF, get the right expressions matrix,
	# the gene annotations, and make the Bioconductor object (ExpressionSet,
	# MAList, or SummarizedExperiment).
	atlasExperimentSummary <- lapply( allAnalytics, function( analytics ) {
		
		# Get the SDRF rows for this analytics object's assays.
		analyticsSDRF <- .createAnalyticsSDRF( analytics, atlasSDRF )

		# Create Bioconductor object, either ExpressionSet, SummarizedExperiment, or MAList.
		biocObject <- .createBiocObject( 
			analytics, 
			analyticsSDRF, 
			atlasExperimentType, 
			experimentAccession, 
			atlasExperimentDirectory 
		)
	})

	atlasExperimentSummary <- SimpleList( atlasExperimentSummary )

	return( atlasExperimentSummary )

}


# .parseSDRF
# 	- Take an SDRF filename and the experiment type from the XML config, and
# 	return a subset of the SDRF containing only the assay names (or ENA runs),
# 	the Characteristics, and FactorValue columns. 
# 	- It tries to include unit columns as well.
# 	- It renames the columns to remove e.g. Characteristics[].
# 	- It removes duplicated columns, e.g. if genotype is a Characteristic and a
# 	Factor.
# 	- It returns the new "SDRF" as a data frame.
.parseSDRF <- function( filename, atlasExperimentType ) {

	# Read in the SDRF file. Set header=FALSE because we don't want the column
	# headings to be made "R-safe" -- this confuses things when we're trying to
	# find the Charactersitic and Factor names. Set stringsAsFactors=FALSE so
	# that we can use grep.
	completeSDRF <- read.delim( filename, header = FALSE, stringsAsFactors = FALSE )

	# Get the Characteristics column indices, and any unit columns next to them.
	charColIndices <- grep( "Characteristics", completeSDRF[ 1, ] )
	charColIndices <- .addUnitCols( charColIndices, completeSDRF )
	
	# Get the Factor column indices, an any unit columns next to them.
	factorColIndices <- grep( "Factor\\s?Value", completeSDRF[ 1, ] )
	factorColIndices <- .addUnitCols( factorColIndices, completeSDRF )

	# Get the index of the Comment[technical replicate group] column, if there
	# is one.
	techRepGroupColIndex <- grep( "Comment\\s?\\[\\s?technical[ _]replicate[ _]group\\s?\\]", completeSDRF[ 1, ] )
	
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

	# If we got a technical replicate group column, add this at the end of the
	# subsetSDRF.
	if( length( techRepGroupColIndex ) > 0 ) {
		subsetSDRF <- cbind( subsetSDRF, completeSDRF[ , techRepGroupColIndex ] )
	}

	# Next thing is to name the columns so they have nice names.
	newColNames <- gsub( "Characteristics\\s?\\[", "", subsetSDRF[1,] )
	newColNames <- gsub( "Factor\\s?Value\\s?\\[", "", newColNames )
	newColNames <- gsub( "Unit\\s?\\[", "", newColNames )
	newColNames <- gsub( "\\s?\\]", "", newColNames )
	newColNames[ 1 ] <- "AssayName"

	# Replace the last column name with one for the technical replicate group
	# column, if there is one.
	if( length( techRepGroupColIndex ) > 0 ) {
		newColNames[ length( newColNames ) ] <- "technical_replicate_group"
	}
	
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
	
	print( subsetSDRF )
	print( newColNames )

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


# .addUnitCols
# 	- Given a vector of column indices and a data frame with the complete SDRF,
# 	return a vector of column indices containing the original ones plus any
# 	Unit[] columns that are next to them.
.addUnitCols <- function( colIndices, SDRF ) {

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


# .createBiocObject
# 	- Takes an Analytics object, experiment type, experiment accession, and
# 	path to directory containing expressions file.
# 	- Returns either an ExpressionSet, MAList, or SummarizedExperiment.
.createBiocObject <- function( analytics, analyticsSDRF, atlasExperimentType, experimentAccession, atlasExperimentDirectory ) {
	
	# Create path to analysis methods file.
	analysisMethodsFile <- paste( experimentAccession, "-analysis-methods.tsv", sep="" )
	analysisMethodsFile <- file.path( atlasExperimentDirectory, analysisMethodsFile )
	
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
		expressions <- read.delim( expressionsFile, header=TRUE, stringsAsFactors=FALSE, row.names = 1 )
		
		# Create an ExpressionSet.
		biocObject <- .createExpressionSet( expressions, analyticsSDRF, analysisMethodsFile )
		
		# Return it.
		return( biocObject )
	}
	# Is this a 2-colour array experiment?
	else if( grepl( "microarray_2colour", atlasExperimentType ) ) {

		# Get the array design.
		arrayDesign <- platform( analytics )
		
		# Create the name of the file with log2 fold-changes.
		mValuesFile <- paste( experimentAccession, "_", arrayDesign, "-log-fold-changes.tsv.undecorated", sep="" )
		mValuesFile <- file.path( atlasExperimentDirectory, mValuesFile )

		# Create the name of the file with average intensities.
		aValuesFile <- paste( experimentAccession, "_", arrayDesign, "-average-intensities.tsv.undecorated", sep="" )
		aValuesFile <- file.path( atlasExperimentDirectory, aValuesFile )
		
		# Read files.
		mValues <- read.delim( mValuesFile, header = TRUE, stringsAsFactors = FALSE, row.names = 1 )
		aValues <- read.delim( aValuesFile, header = TRUE, stringsAsFactors = FALSE, row.names = 1 )

		# Create an MAList
		maList <- list(
			genes = rownames( mValues ),
			M = mValues,
			A = aValues
		)

		biocObject <- new( "MAList", maList )

		return( biocObject )
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
		expressions <- read.delim( expressionsFile, header=TRUE, stringsAsFactors=FALSE, row.names = 1 )
		
		# Create a SummarizedExperiment.
		biocObject <- .createSummarizedExperiment( expressions, analyticsSDRF, analysisMethodsFile )
	}
	# Otherwise, don't recognise this experiment type so die.
	else {
		stop( paste( "Don't know how handle experiment type \"", atlasExperimentType, "\". Cannot continue.", sep="" ) )
	}
	
	# Return the data frame with expressions.
	return( biocObject )
}


# .createAnalyticsSDRF
# 	- Take an Analytics object and a parsed SDRF.
# 	- Return a new SDRF data frame with just the assay_names from the Analytics
# 	object, as well as a column denoting which assay group each assay belongs
# 	to.
.createAnalyticsSDRF <- function( analytics, atlasSDRF ) {

	# Get the assay groups.
	assayGroups <- assay_groups( analytics )

	# Get the SDRF rows for these assay groups.
	assayGroupSDRFs <- lapply( assayGroups, function( assayGroup ) {
		
		# Get the assay names.
		assayNames <- assay_names( assayGroup )

		# Sanity checking.
		if( length( assayNames ) == 0 ) { 
			stop( "ERROR - Did not find any assay names for an assay group. Cannot continue." )
		}

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


# .createExpressionSet
# 	- Take a data frame of normalized expressions and a parsed SDRF.
# 	- Return an ExpressionSet object.
.createExpressionSet <- function( expressions, analyticsSDRF, analysisMethodsFile ) {

	# Only select columns with assay names in our SDRF -- these are the
	# ones that passed QC and are still in the XML.
	expressions <- expressions[ , rownames( analyticsSDRF ) ]

	# Turn data frame into matrix.
	expressionsMatrix <- as.matrix( expressions )
	
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
	analysisMethodsList <- .readArrayAnalysisMethods( analysisMethodsFile )
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


.createSummarizedExperiment <- function( expressions, analyticsSDRF, analysisMethodsFile ) {
	
	# Only select columns with assay names in our SDRF.
	expressions <- expressions[ , rownames( analyticsSDRF ) ]

	# Turn data frame into matrix.
	expressionsMatrix <- as.matrix( expressions )

	# Turn SDRF into a DataFrame.
	analyticsSDRF <- DataFrame( analyticsSDRF )

	# Get analysis methods
	analysisMethodsList <- .readSeqAnalysisMethods( analysisMethodsFile )

	# Create SummarizedExperiment
	summarizedExperiment <- SummarizedExperiment( assays = SimpleList( counts = expressionsMatrix ), colData = analyticsSDRF, exptData = analysisMethodsList )

	# Return it.
	return( summarizedExperiment )
}


.readArrayAnalysisMethods <- function( analysisMethodsFile ) {
	
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


.readSeqAnalysisMethods <- function( analysisMethodsFile ) {
	
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
