#!/usr/bin/env Rscript

# diffAtlas_DE_limma.R
# Microarray differential expression statistics computation for the
# Differential Atlas.

suppressMessages( library( limma ) )
suppressMessages( library( Biobase ) )
suppressMessages( library( genefilter ) )

suppressMessages( library( ExpressionAtlas ) )

# diffAtlas_DE_limma()
# - Differential expression analysis (2-group comparison) using limma.
# Arguments:
# 	- normExprsFile <- matrix of normalized and summarized expression values.
#	- refAssays <- comma-separated list of assay accessions in reference assay group.
#	- testAssays <- comma-separated list of assay accessions in test assay group.
#	- resFile <- filename for results.
#	- plotDataFile <- filename for data for MvA plot.
#	- aValuesFile <- filename for matrix of A-values; only needed for 2-colour designs. Default is NULL.
#diffAtlas_DE_limma <<- function( normExprsFile, refAssays, testAssays, resFile, plotDataFile, aValuesFile = NULL ) {
diffAtlas_DE_limma <<- function( expAcc, xmlConfigFileName, atlasProcessingDirectory ) {

	e <- try({

		# First we need to parse the config file.
		cat( paste( "Reading XML config from", xmlConfigFilename, "..." ) )
		experimentConfig <- parseAtlasConfig( xmlConfigFilename )
		cat( paste( "Successfully read XML config." ) )

		# Get the experiment type.
		experimentType <- experimentConfig$experimentType

		# Check that this is not an RNA-seq experiment.
		if( !grepl( "array", experimentType ) ) {
			stop( paste(
						"Experiment type",
						experimentType,
						"does not look like a microarray experiment. Cannot continue" 
			) )
		}

		# Get the list of analytics objects from the config.
		allAnalytics <- experimentConfig$allAnalytics
			
		# Steps are different for 1-colour and 2-colour data.
		if( grepl( "1colour", experimentType ) ) {

			run_one_colour_analysis( expAcc, allAnalytics, atlasProcessingDirectory )

		} else if( grepl( "2colour", experimentType ) ) {

			run_two_colour_analysis( expAcc, allAnalytics, atlasProcessingDirectory )

		} else {
			cat( paste(
					   "Experiment type",
					   experimentType,
					   "not recognised. Cannot continue."
			) )
		}

	} )
}
		


select_assay_groups_for_contrast <- function( expAssayGroups, refAssayGroupID, testAssayGroupID ) {

	# Get the assay groups.
	refAssayGroup <- expAssayGroups[[ refAssayGroupID ]]
	testAssayGroup <- expAssayGroups[[ testAssayGroupID ]]

	# Return them in a list.
	return( 
		   list( 
				 refAssayGroup = refAssayGroup,
				 testAssayGroup = testAssayGroup
				 ) 
	)
}


make_biorep_annotations <- function( contrastAssayGroups, expBatchEffects, twoColour ) {

	if( missing( twoColour ) ) { twoColour <- 0 }

	refAssaysToBioReps <- make_assays_to_bioreps_df( biological_replicates( contrastAssayGroups$refAssayGroup ), twoColour )
	testAssaysToBioReps <- make_assays_to_bioreps_df( biological_replicates( contrastAssayGroups$testAssayGroup ), twoColour )
	
	# Mappings of assay names to biological replicate names. Biological
	# replicate names are most commonly the same as the assay names except in
	# the case of technical replicates, in which case the biological replicate
	# name is the technical replicate group ID.
	allAssaysToBioReps <- rbind( refAssaysToBioReps, testAssaysToBioReps )

	# Remove the row names as they're not useful.
	rownames( allAssaysToBioReps ) <- NULL

	# Get the contrast assay names to subset the batch effect lists with.
	contrastAssayNames <- allAssaysToBioReps$AssayName
	
	# Make the column for the group (i.e. test or reference) that the assays
	# belong to.
	groupsCol <- c( rep( "ref", nrow( refAssaysToBioReps ) ), rep( "test", nrow( testAssaysToBioReps ) ) )

	# Add the groups column to the data frame.
	allAssaysToBioReps$groups <- groupsCol

	if( length( expBatchEffects ) > 0 ) {
	
		# Get a list of data frames mapping assay names to batch effects.
		allBatchEffectDfs <- lapply( expBatchEffects, function( batchEffect ) {
			batch_effect_to_df( batchEffect, contrastAssayNames )
		})
		
		# Combine all batch effect data frames into one.
		batchEffectsDf <- do.call( "cbind", allBatchEffectDfs )

		# Sort the rows of batchEffectsDf in the same order as the assay names in
		# allAssaysToBioReps.
		batchEffectsDf <- batchEffectsDf[ allAssaysToBioReps$AssayName , ]

		# Combine the batch effects data frame and the assays-to-bioreps data frame
		# to make the final data frame.
		return( cbind( allAssaysToBioReps, batchEffectsDf ) )
	
	} else {
		
		return( allAssaysToBioReps )
	}
}


make_assays_to_bioreps_df <- function( bioReps, twoColour ) {

	if( missing( twoColour ) ) { twoColour <- 0 }

	assaysToBioRepsList <- lapply( bioReps, function( bioRep ) {

		assayNames <- biorep_assay_names( bioRep )
		techRepId <- technical_replicate_id( bioRep )

		if( length( techRepId ) > 0 ) {
			
			if( twoColour ) {
				
				dyeNames <- sub( ".*(Cy\\d)$", "\\1", assayNames )
				
				if( length( unique( dyeNames ) ) > 1 ) {
					stop( "Technical replicates with different dye names are not allowed." )
				}

				dyeName <- unique( dyeNames )

				techRepId <- paste( techRepId, dyeName, sep = "." )
			}

			data.frame( AssayName = assayNames, BioRepName = techRepId, stringsAsFactors = FALSE )

		} else {
			
			data.frame( AssayName = assayNames, BioRepName = assayNames, stringsAsFactors = FALSE )
		}
	})

	assaysToBioRepsDf <- do.call( "rbind", assaysToBioRepsList )

	return( assaysToBioRepsDf )
}


batch_effect_to_df <- function( batchEffect, contrastAssayNames ) {

	assayBatches <- batches( batchEffect )
	effectName <- effect_name( batchEffect )

	# Get a list of data frames with the assays from each batch.
	batchDataFrames <- lapply( names( assayBatches ), function( batchName ) {
		
		assayNames <- assayBatches[[ batchName ]]

		data.frame( AssayNames = assayNames, batch = batchName, stringsAsFactors = FALSE )
	} )

	# Merge the data frames in to one data frame.
	allBatches <- do.call( "rbind", batchDataFrames )

	# Change the column heading from "batch" to the actual effect name.
	colnames( allBatches ) <- c( "AssayName", effectName )

	# Subset only the contrast's assay names.
	allBatches <- allBatches[ which( allBatches$AssayName %in% contrastAssayNames ), ]

	# Re-order the rows based on the assay names.
	allBatches <- allBatches[ order( allBatches$AssayName ), ]

	rownames( allBatches ) <- allBatches$AssayName
	allBatches$AssayName <- NULL

	return( allBatches )
}


add_tech_rep_averages <- function( expressionData, bioRepAnnotations, techRepGroupIDs, twoColour ) {

	if( missing( twoColour ) ) { twoColour <- 0 }

	# If we have 2-colour data, we only want to average between the assays
	# once, but the assay names are in the data frame twice, once for each dye.
	if( twoColour ) {
		
		assayNamesNoCy <- gsub( ".Cy\\d$", "", bioRepAnnotations$AssayName )
		bioRepNamesNoCy <- gsub( ".Cy\\d$", "", bioRepAnnotations$BioRepName )
		
		uniqueTechRepGroupIDs <- unique( gsub( ".Cy\\d$", "", techRepGroupIDs ) )

		techRepAveragesList <- lapply( uniqueTechRepGroupIDs, function( techRepGroupID ) {

			techRepIndices <- which( bioRepNamesNoCy == techRepGroupID )

			techRepAssayNames <- assayNamesNoCy[ techRepIndices ]

			uniqueTechRepAssayNames <- unique( techRepAssayNames )

			techRepColumns <- expressionData[ , uniqueTechRepAssayNames ]

			techRepAverages <- apply( techRepColumns, 1, function( x ) { mean( x ) } )

			techRepAveragesDF <- structure( data.frame( techRepAverages ), names = techRepGroupID )
		} )
	}
    else {

		# Create a list with one element per column of averages.
		techRepAveragesList <- lapply( techRepGroupIDs, function( techRepGroupID ) {

			# Get the row indices for the technical replicates.
			techRepRows <- which( bioRepAnnotations$BioRepName == techRepGroupID )

			# Get the assay names on these rows.
			techRepAssayNames <- bioRepAnnotations$AssayName[ techRepRows ]

			# Get the columns of normalized data for these assays.
			techRepColumns <- expressionData[ , techRepAssayNames ]

			# Get the mean of each row.
			techRepAverages <- apply( techRepColumns, 1, function( x ) { mean( x ) } )

			techRepAveragesDF <- structure( data.frame( techRepAverages ), names = techRepGroupID )
		} )
	}

	# Turn the list into a data frame.
	techRepAveragesDF <- do.call( "cbind", techRepAveragesList )

	# Next need to remove the columns for the technical replicate assays from
	# the normalized data frame.
	# First get all the assay names we want to remove.
	allTechRepRows <- which( bioRepAnnotations$BioRepName %in% techRepGroupIDs )

	# Now get the assay names for those rows.
	if( twoColour ) {
		allTechRepAssays <- unique( assayNamesNoCy[ allTechRepRows ] )

	} else { 
		allTechRepAssays <- bioRepAnnotations$AssayName[ allTechRepRows ]
	}

	# Remove the columns for these assay names from the normalized data frame.
	expressionData <- expressionData[ , !( colnames( expressionData ) %in% allTechRepAssays ), drop = FALSE ]

	# Add the columns of averages instead.
	expressionData <- cbind( expressionData, techRepAveragesDF )

	return( expressionData )
}
	

filter_and_adjust_pvalues <- function( rawPvalues, normDataRowVars ) {

	# Independent filtering.
	# Make a data frame containing the row variances of the normalized data
	# matrix, and the unadjusted p-values.
	filterData <- data.frame( rowvar = normDataRowVars, test = rawPvalues )
	
	# theta is a vector of numbers from 0 to 0.8 in increments of 0.02.
	# These values represent the proportion of the lower end of the data to
	# filter out. We will use them to find out how many true null
	# hypotheses we will be rejecting, if we filter out the proportion
	# corresponding to each of them.
	theta = seq( from=0, to=1.0, by=0.02 )
	
	# Work out adjusted p-values after filtering out each proportion of
	# data specified in theta.
	filteredAdjustedPvalues <- filtered_p( 
		filter = filterData$rowvar,	# use the row variances as the filter statistic.
		test = filterData$test,	# the unadjusted p-values.
		theta = theta,	# the range of filtering proportions.
		method = "BH"	# use Benjamini-Hochberg for p-value adjustment.
	)

	# filteredAdjustedPvalues is a matrix of adjusted p-values, with a
	# column for each proportion from theta.  For each column, count how
	# many rejections of the null hypothesis we will make, using the FDR
	# threshold of 0.05.
	numRej <- colSums( filteredAdjustedPvalues < 0.05, na.rm=TRUE )
	
	# Find the index of the column that had the most rejections of the null
	# hypothesis.
	maxRejectionsIndex <- which.max( numRej )

	# Return the column of adjusted p-values at this index to fitEbayes.
	return( filteredAdjustedPvalues[ , maxRejectionsIndex ] )
}


make_eset_for_contrast <- function( normalizedData, bioRepAnnotations ) {

	# Turn normalized data into a matrix.
	exprsForContrast <- as.matrix( normalizedData )
	
	# Create a new AssayData object containing the normalized data.
	expressionData <- assayDataNew( storage.mode = "lockedEnvironment", exprs = exprsForContrast )
	
	# Add the featureNames and sampleNames attributes.
	featureNames( expressionData ) <- rownames( exprsForContrast )
	sampleNames( expressionData ) <- colnames( exprsForContrast )
	
	# Create the sample annotations AnnotatedDataFrame object.
	sampleAnnotations <- new( "AnnotatedDataFrame", data = bioRepAnnotations )

	# Create an AnnotatedDataFrame object with the feature data (design element
	# names).
	featureData <- new( "AnnotatedDataFrame", data = data.frame( designElements = rownames( normalizedData ) ) )
	featureNames( featureData ) <- rownames( normalizedData )

	return( new(
				"ExpressionSet",
				assayData = expressionData,
				phenoData = sampleAnnotations,
				featureData = featureData
	) )
}


check_twocolour_assaynames <- function( contrastAssayGroups ) {

	# Get the test and reference assay groups.
	testAssayGroup <- contrastAssayGroups$testAssayGroup
	refAssayGroup <- contrastAssayGroups$refAssayGroup
	
	# Get the assay names from each assay group and remove the dye names.
	testAssaysNoCy <- gsub( ".Cy\\d$", "", assay_names( testAssayGroup ) )
	refAssaysNoCy <- gsub( ".Cy\\d$", "", assay_names( refAssayGroup ) )
	
	# Sort the assay names so they're in the same order.
	testAssaysNoCy  <- sort( testAssaysNoCy )
	refAssaysNoCy <- sort( refAssaysNoCy )

	# Die if the test and reference assay names without dye names are not identical.
	if( any( testAssaysNoCy != testAssaysNoCy ) ) {
		stop(  "Differing assay names found in test and reference assay groups after removing \".Cy3\" and \".Cy5\". Please verify this experiment has a two-colour design." )
	}
}


select_twocolour_columns <- function( allData, contrastAssayGroups ) {

	assayNamesNoCy <- gsub( ".Cy\\d$", "", assay_names( contrastAssayGroups$testAssayGroup ) )
	
	wantedCols <- allData[ , assayNamesNoCy ]
	
	return( wantedCols )
}


run_one_colour_analysis <- function( expAcc, allAnalytics, atlasProcessingDirectory ) {

	# Go through the analytics and do the analysis...
	invisible( sapply( allAnalytics, function( analytics ) {

		# Get the platform (array design).
		arrayDesign <- platform( analytics )

		# Create the normalized data file name.
		normalizedDataFilename <- paste( 
										expAcc, 
										"_", 
										arrayDesign, 
										"-normalized-expressions.tsv.undecorated", 
										sep = "" 
										)

		normalizedDataFilename <- file.path( atlasProcessingDirectory, normalizedDataFilename )

		# Read in the normalized data.
		normalizedData <- read.delim( 
									 normalizedDataFilename, 
									 header = TRUE, 
									 stringsAsFactors = FALSE,
									 row.names = 1
									 )

		# Get the contrasts, assay groups, and batch effects.
		expContrasts <- atlas_contrasts( analytics )
		expAssayGroups <- assay_groups( analytics )
		expBatchEffects <- batch_effects( analytics )

		# Go through the contrasts...
		sapply( expContrasts, function( expContrast ) {
			
			# Get the contrast info.
			contrastID <- contrast_id( expContrast )
			contrastName <- contrast_name( expContrast )
			refAssayGroupID <- reference_assay_group_id( expContrast )
			testAssayGroupID <- test_assay_group_id( expContrast )

			# Get the two assay groups for this contrast as a list.
			contrastAssayGroups <- select_assay_groups_for_contrast( 
																 expAssayGroups, 
																 refAssayGroupID, 
																 testAssayGroupID 
																 )
			
			# Create the data frame containing annotations for each
			# biological replicate, grouping technical replicates.
			bioRepAnnotations <- make_biorep_annotations( contrastAssayGroups, expBatchEffects )

			# Next, subset the normalized data frame so that it
			# contains only the columns of data for this contrast.
			normalizedData <- normalizedData[ , bioRepAnnotations$AssayName ]
			
			# Get the technical replicate group IDs -- these are the
			# ones in the BioRepName column that are duplicated.
			techRepGroupIDs <- unique( bioRepAnnotations$BioRepName[ which( duplicated( bioRepAnnotations$BioRepName ) ) ] )

			# If there are any technical replicates...
			if( length( techRepGroupIDs ) > 0 ) {
				# Replace the columns for technical replicates with their averages in the normalized data.
				normalizedData <- add_tech_rep_averages_1colour( normalizedData, bioRepAnnotations, techRepGroupIDs )
			}
			
			# Remove the AssayName column from the annotations data frame now we no longer need it.
			bioRepAnnotations$AssayName <- NULL

			# Remove duplicated technical replicate rows, if any.
			duplicatedRowIndices <- which( duplicated( bioRepAnnotations$BioRepName ) )
			if( length( duplicatedRowIndices ) > 0 ) {
				bioRepAnnotations <- bioRepAnnotations[ -duplicatedRowIndices, ]
			}

			# Add the biological replicate names as the row names.
			rownames( bioRepAnnotations ) <- bioRepAnnotations$BioRepName

			# Make the column names R-safe.
			colnames( bioRepAnnotations ) <- make.names( colnames( bioRepAnnotations ) )

			# Put the columns of the normalized data into the same
			# order as the rows in the bio rep annotations.
			normalizedData <- normalizedData[ , rownames( bioRepAnnotations ) ]
			
			# Create the ExpressionSet object.
			esetForContrast <- make_eset_for_contrast( normalizedData, bioRepAnnotations )

			# Make the formula to create the design matrix, adding
			# batch effects if necessary.
			# If we have batch effects...
			if( length( expBatchEffects ) > 0 ) {
				
				# Check that some batch effects exist in the
				# ExpressionSet's info, if not something went wrong.
				if( ncol( pData( esetForContrast ) ) < 3 ) {
					stop( "Expecting batch effect information in ExpressionSet but none was found." )
				}
				
				# Get the batch effect names -- these should be the same as
				# the column headings in the ExpressionSet's info.
				batchEffectNames <- make.names( sapply( expBatchEffects, function( batchEffect ) { effect_name( batchEffect ) } ) )

				# Check that the batch effect names match the column
				# headings in the ExpressionSet.
				if( !all( ( batchEffectNames %in% colnames( pData( esetForContrast ) ) ) ) ) {
					nameString <- paste( batchEffectNames, collapse = ", " )
					stop( paste( "Did not find all batch effect names in ExpressionSet object. Names are:", nameString ) )
				}

				# Create the formula string.
				formulaString <- paste( "~", paste( batchEffectNames, collapse = " + " ), "+ groups" )
			
			} else {

				# If we don't have any batch effects, create the formula string without them.
				formulaString <- "~ groups"
			}

			# Now create the design matrix.
			designMatrix <- model.matrix( as.formula( formulaString ) , data = esetForContrast )

			# Fit the linear model.
			fit <- lmFit( esetForContrast, design = designMatrix )
			
			# Run empirical Bayes method to get differential expression statistics.
			fit <- eBayes( fit )
			
			# Check that we have a column called "groupstest" in the
			# results. If not then things didn't work out as expected.
			stopifnot( "groupstest" %in% colnames( fit$p.value ) )
			
			# Adjust the p-values and perform independent filtering, add to the fit object.
			fit$adjPvals <- filter_and_adjust_pvalues( rowVars( normalizedData ), fit$p.value[ , "groupstest" ] )

			# results for heatmap matrix display:
			contrastResults <- data.frame(
				designElements = featureNames( esetForContrast ), 
				adjPval = fit$adjPvals, 
				t = fit$t[ , "groupstest" ], 
				logFC = fitEbayes$coefficients[ , "groupstest" ]
			)

			# results to be used for MvA plot creation (as above but with average intensities and without t-stats):
			plotData <- data.frame(
				designElements = featureNames(esetForContrast), 
				adjPval = fit$adjPvals, 
				logFC = fitEbayes$coefficients[ , "groupstest" ], 
				avgExpr = fitEbayes$Amean[ , "groupstest" ]
			)
			
			# Create filenames to write to.
			# Stats results:
			resFile <- paste( expAcc, contrastID, "analytics", "tsv", sep = "." )
			resFile <- file.path( Sys.getenv( "HOME" ), "tmp", resFile )
			# Data for MvA plot:
			plotDataFile <- paste( expAcc, contrastID, "plotdata", "tsv", sep = "." )
			plotDataFile <- file.path( Sys.getenv( "HOME" ), "tmp", plotDataFile )
			
			# Write the files.
			write.table( contrastResults, file=resFile, row.names=FALSE, quote=FALSE, sep="\t" )
			write.table( plotData, file=plotDataFile, row.names=FALSE, quote=FALSE, sep="\t" )

		} )

	} ) )
}


run_two_colour_analysis <- function( expAcc, allAnalytics, atlasProcessingDirectory ) {

	# Go through the analytics and do the analysis...
	invisible( sapply( allAnalytics, function( analytics ) {

		# Get the platform (array design).
		arrayDesign <- platform( analytics )

		# Create the log2 fold-changes file name.
		logfcsFilename <- paste( 
								expAcc, 
								"_", 
								arrayDesign, 
								"-log-fold-changes.tsv.undecorated", 
								sep = "" 
								)

		logfcsFilename <- file.path( atlasProcessingDirectory, logfcsFilename )

		# Read in the log fold-changes.
		logFoldChanges <- read.delim( 
								 logfcsFilename, 
								 header = TRUE, 
								 stringsAsFactors = FALSE,
								 row.names = 1
								 )

		# Create the average intensities file name.
		aveIntsFilename <- paste(
								 expAcc,
								 "_",
								 arrayDesign,
								 "-average-intensities.tsv.undecorated",
								 sep = ""
								 )
		aveIntsFilename <- file.path( atlasProcessingDirectory, aveIntsFilename )

		# Read in the average intensities.
		averageIntensities <- read.delim(
										 aveIntsFilename,
										 header = TRUE,
										 stringsAsFactors = FALSE,
										 row.names = 1
										 )

		# Get the contrasts, assay groups, and batch effects.
		expContrasts <- atlas_contrasts( analytics )
		expAssayGroups <- assay_groups( analytics )
		expBatchEffects <- batch_effects( analytics )
		# We shouldn't have any batch effect info for 2-colour array data as
		# can't handle this yet.
		if( length( expBatchEffects ) > 0 ) {
			stop( "Cannot handle batch effects for 2-colour microarray data." )
		}

		# Go through the contrasts...
		sapply( expContrasts, function( expContrast ) {
			
			# Get the contrast info.
			contrastID <- contrast_id( expContrast )
			contrastName <- contrast_name( expContrast )
			refAssayGroupID <- reference_assay_group_id( expContrast )
			testAssayGroupID <- test_assay_group_id( expContrast )

			# Get the two assay groups for this contrast as a list.
			contrastAssayGroups <- select_assay_groups_for_contrast( 
																 expAssayGroups, 
																 refAssayGroupID, 
																 testAssayGroupID 
																 )

			# Make sure the assay names for test and ref groups are identical
			# once the dye names are removed. Dies here if not.
			check_twocolour_assaynames( contrastAssayGroups )

			# Select only the columns we need from the log fold-changes and average intensities.
			logFoldChanges <- select_twocolour_columns( logFoldChanges, contrastAssayGroups )
			averageIntensities <- select_twocolour_columns( averageIntensities, contrastAssayGroups )
			
			# Make the biological replicate annotations data frame.
			bioRepAnnotations <- make_biorep_annotations( contrastAssayGroups, expBatchEffects, 1 )
		
			# Get the technical replicate group IDs -- these are the
			# ones in the BioRepName column that are duplicated.
			techRepGroupIDs <- unique( bioRepAnnotations$BioRepName[ which( duplicated( bioRepAnnotations$BioRepName ) ) ] )

			if( length( techRepGroupIDs ) > 0 ) {
				# Replace the columns for technical replicates with their
				# averages in the log fold-changes and average intensities.
				logFoldChanges <- add_tech_rep_averages( logFoldChanges, bioRepAnnotations, techRepGroupIDs, 1 )
				averageIntensities <- add_tech_rep_averages( averageIntensities, bioRepAnnotations, techRepGroupIDs, 1 )
			}

			# Remove the AssayName column from the annotations data frame now we no longer need it.
			bioRepAnnotations$AssayName <- NULL

			# Remove duplicated technical replicate rows, if any.
			duplicatedRowIndices <- which( duplicated( bioRepAnnotations$BioRepName ) )
			if( length( duplicatedRowIndices ) > 0 ) {
				bioRepAnnotations <- bioRepAnnotations[ -duplicatedRowIndices, ]
			}

			# Make the column names R-safe.
			colnames( bioRepAnnotations ) <- make.names( colnames( bioRepAnnotations ) )

			cy3indices <- grep( "Cy3", bioRepAnnotations$BioRepName )
			cy5indices <- grep( "Cy5", bioRepAnnotations$BioRepName )

			cy3bioreps <- sub( ".Cy\\d$", "", bioRepAnnotations$BioRepName )[ cy3indices ]
			cy5bioreps <- sub( ".Cy\\d$", "", bioRepAnnotations$BioRepName )[ cy5indices ]

			cy3groups <- bioRepAnnotations$groups[ cy3indices ]
			cy5groups <- bioRepAnnotations$groups[ cy5indices ]

			cy3df <- data.frame( BioRepName = cy3bioreps, groups = cy3groups )
			cy5df <- data.frame( BioRepName = cy5bioreps, groups = cy5groups )

			cy5df <- cy5df[ order( cy5df$BioRepName ), ]
			cy3df <- cy3df[ order( cy3df$BioRepName ), ]

			stopifnot( all( cy3df$BioRepName == cy5df$BioRepName ) )

			targetsDF <- data.frame( Cy3 = cy3df$groups, Cy5 = cy5df$groups )
			rownames( targetsDF ) <- cy3df$BioRepName

			designMatrix <- modelMatrix( targets, ref = "ref" )

			logFCsForContrast <- logFoldChanges[ , rownames( targetsDF ) ]
			avgIntsForContrast <- averageIntensities[ , rownames( targetsDF ) ]

			maList <- list(
						   genes = rownaes( logFCsForContrast ),
						   M = logFCsForContrast,
						   A = avgIntsForContrast
						   )
			maList <- new( "MAList", maList )

			fit <- lmFit( maList, designMatrix )

			fit <- eBayes( fit )

			# Adjust the p-values and perform independent filtering, add to the fit object.
			fit$adjPvals <- filter_and_adjust_pvalues( rowVars( logFCsForContrast ), fit$p.value )

			contrastResults <- data.frame(
										  designElements = fit$genes,
										  adjPval = fit$adjPvals,
										  t = fit$t[ , 1 ],
										  logFC = fit$coefficients[ , 1 ]
										  )
			plotData <- data.frame(
								   designElements = fit$genes,
								   adjPval = fit$adjPvals,
								   logFC = fit$coefficients[ , 1 ],
								   avgExpr = fit$Amean
								   )

			# Create filenames to write to.
			# Stats results:
			resFile <- paste( expAcc, contrastID, "analytics", "tsv", sep = "." )
			resFile <- file.path( Sys.getenv( "HOME" ), "tmp", resFile )
			# Data for MvA plot:
			plotDataFile <- paste( expAcc, contrastID, "plotdata", "tsv", sep = "." )
			plotDataFile <- file.path( Sys.getenv( "HOME" ), "tmp", plotDataFile )
			
			# Write the files.
			write.table( contrastResults, file=resFile, row.names=FALSE, quote=FALSE, sep="\t" )
			write.table( plotData, file=plotDataFile, row.names=FALSE, quote=FALSE, sep="\t" )
		
		} )

	} ) )
}

		


	
###################
# RUN MAIN FUNCTION
###################

# Run with arguments if there are any, otherwise don't do anything. Having this
# lets us source this file in an R session and load all the functions so we can
# run them separately if desired.
args <- commandArgs( TRUE )
if( length( args ) > 0) {
	do.call( diffAtlas_DE_limma, as.list( args ) )
}
