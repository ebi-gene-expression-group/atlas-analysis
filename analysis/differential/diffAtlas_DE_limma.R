#!/usr/bin/env Rscript

# diffAtlas_DE_limma.R
# Microarray differential expression statistics computation for the
# Differential Atlas.

library( limma )
library( Biobase )
library( genefilter )

# diffAtlas_DE_limma()
# - Differential expression analysis (2-group comparison) using limma.
# Arguments:
# 	- normExprsFile <- matrix of normalized and summarized expression values.
#	- refAssays <- comma-separated list of assay accessions in reference assay group.
#	- testAssays <- comma-separated list of assay accessions in test assay group.
#	- resFile <- filename for results.
#	- plotDataFile <- filename for data for MvA plot.
#	- aValuesFile <- filename for matrix of A-values; only needed for 2-colour designs. Default is NULL.
diffAtlas_DE_limma <<- function( normExprsFile, refAssays, testAssays, resFile, plotDataFile, aValuesFile = NULL ) {

	e <- try({

		twoColour <- NULL
		# Set twoColour if there is an aValueFile.
		if( !is.null( aValuesFile ) ) twoColour <- 1
		
		# Read expressions (M-values/log-fold-changes for 2-colour).
		print( paste( "Reading", normExprsFile ) )
		# R's reserved characters are changed to "." because check.names=TRUE
		# in the read.delim() function.
		normExprs <- read.delim( normExprsFile, stringsAsFactors = FALSE )
		rownames( normExprs ) <- normExprs[ , 1 ]
		normExprs[ , 1 ] <- NULL

		# Make vectors of assay accessions.
		# Assay groups come in like e.g.:
		# 	assay 1<SEP>assay 2<SEP>assay 3
		# Split on "<SEP>" -- we use this as a separator because it is very
		# unlikely to crop up in actual assay names.
		refAssays <- unlist( strsplit( refAssays, "<SEP>" ) )
		testAssays <- unlist( strsplit( testAssays, "<SEP>" ) )
		

		# Assays that are technical replicates are separated by e.g.
		# "<T1_SEP>", where the number after "T" is the number from the technical
		# replicate group comment in the SDRF.
		
		# Get vectors of the technical replicate group(s), if any.
		refAssays_techreps <- refAssays[ grep( "<T\\d+_SEP>", refAssays ) ]
		testAssays_techreps <- testAssays[ grep( "<T\\d+_SEP>", testAssays ) ]
		
		# If there are some technical replicates...

		# Pass normExprs and tech rep assay names to function. This will
		# replace columns for the tech rep assays in normExprs with their
		# mean values, and then return normExprs. Do this for the reference
		# assays and/or the test assays. Also replace assay names of tech
		# rep assays with tech rep group names. We pass the twoColour
		# variable to the functions used here as well, because for
		# two-colour array data we need to add back the "Cy3" or "Cy5" to
		# the technical replicate group names after we've created them.
		if( length( refAssays_techreps ) > 0 ) {
			
			# Replace columns for tech rep assays with their averages.
			normExprs <- addTechRepAverages( normExprs, refAssays_techreps, twoColour )
			
			# Replace names of tech rep assays with tech rep group IDs to
			# match columns added above.
			refAssays <- replaceTechRepAssayNames( refAssays, twoColour )
		}
	
		# If this is a one-colour design, also look at the test assays.
		if( is.null( twoColour ) ) {
		
			if( length( testAssays_techreps ) > 0 ) {
				
				# Replace columns for tech rep assays with their averages.
				normExprs <- addTechRepAverages( normExprs, testAssays_techreps, twoColour )
				
				# Replace names of tech rep assays with tech rep group IDs to
				# match columns added above.
				testAssays <- replaceTechRepAssayNames( testAssays, twoColour )
			}
		} 
		# If this is two-colour data, still need to replace the test assay
		# names, but don't need to average columns as that's already done.
		else {
			
			if( length( testAssays_techreps ) > 0 ) {

				# Replace names of tech rep assays with tech rep group IDs to
				# match columns added above.
				testAssays <- replaceTechRepAssayNames( testAssays, twoColour )
			}
		}
			
		# Use make.names to convert reserved characters to "." so that they
		# match the column headings from reading in the normalized expressions
		# above. 
		refAssays <- make.names( refAssays )
		testAssays <- make.names( testAssays )

		# Provisions for 2-colour designs:
		#	- If aValuesFile is not NULL, this means it's a 2-colour design.
		#	- Create an MAList object out of the M-values in normExprs and
		#	A-values in aValues. This is essentially equivalent to the
		#	ExpressionSet created for Affymetrix data.
		#	- Fit LM, get DE stats. Contrast is implied, do not need to specify contrasts matrix.
		#	- Write results and data for MvA plot.
		if( !is.null( twoColour ) ) {

			# Sort ref and test assays so they are in the same order.
			refAssays <- sort( refAssays )
			testAssays <- sort( testAssays )
			
			# Read in A-values files (average intensities).
			aValues <- read.delim( aValuesFile, stringsAsFactors = FALSE )
			rownames( aValues ) <- aValues[ , 1 ]
			aValues[ , 1 ] <- NULL

			# If there are technical replicates, need to go through the same
			# process as above, replacing tech rep assay columns with their
			# averages, for the reference assays and the test assays.
			if( length( refAssays_techreps > 0 ) ) {
				
				# Replace columns in aValues for tech rep assays with their averages.
				aValues <- addTechRepAverages( aValues, refAssays_techreps, 1 )
			}

			# Strip off .Cy* endings to check the assay names are the same and
			# to subset columns from normExprs and aValues.
			refAssaysNoCy <- gsub( ".Cy\\d$", "", refAssays )
			testAssaysNoCy <- gsub( ".Cy\\d$", "", testAssays )

			# [double-]check that we have the same assay names for test and ref
			# without .Cy* endings.
			# Compare sorted values at each index position using
			# any(). If there are any mismatches this returns TRUE.
			if( any( refAssaysNoCy != testAssaysNoCy ) ) {
				stop( "Differing assay names found in test and reference assay groups after removing \".Cy3\" and \".Cy5\". Please verify this experiment has a two-colour design." )
			}
			
			# Create "targets" dataframe for design matrix.
			# Vector of ref and test for Cy3-labelled samples.
			cy3 <- rep( NA, length( refAssays ) )
			cy3[ grep( ".Cy3$", refAssays ) ] <- "ref"
			cy3[ grep( ".Cy3$", testAssays ) ] <- "test"
			# Vector of ref and test for Cy5-labelled samples.
			cy5 <- rep( NA, length( refAssays ) )
			cy5[ grep( ".Cy5$", refAssays ) ] <- "ref"
			cy5[ grep( ".Cy5$", testAssays ) ] <- "test"
			# Check in case any of cy3 or cy5 are still NA -- something's wrong if so.
			if( any( is.na( cy3 ) ) || any( is.na( cy5 ) ) ) {
				stop( "Some assays do not have an accession ending with \".Cy3\" or \".Cy5\". Please verify this experiment has a two-colour design." )
			}
			
			# Create "targets" dataframe like:
			#	-----------
			#		Cy3	Cy5
			#	A1	ref	test
			#	A2	test	ref
			#	A3	ref	test
			#	A4	test	ref
			#	-----------
			targets <- data.frame( Cy3 = cy3, Cy5 = cy5 )
			
			# Rownames are assay names and should be the same as the column
			# headings in dataframes of M-values and A-values.
			rownames( targets ) <- refAssaysNoCy

			# Design matrix
			designMatrix <- modelMatrix( targets, ref="ref" )
			
			# Subset M-values and A-values
			mValuesForContrast <- normExprs[ , refAssaysNoCy ]
			aValuesForContrast <- aValues[ , refAssaysNoCy ]

			# Verify that rownames (i.e. design element IDs) are the same (and
			# in the same order) for subsetted M-values and A-values.
			if( any( rownames( mValuesForContrast ) != rownames( aValuesForContrast ) ) ) {
				stop( "Different design element IDs found in M-values and A-values. Please verify they are from the same experiment." )
			}

			# Create MAList object.
			# First need to make a normal list with "genes" (design element
			# IDs), "M", and "A" elements.
			maList <- list(
				genes = rownames( mValuesForContrast ), 
				M = mValuesForContrast, 
				A = aValuesForContrast
			)
			# Make it into an "MAList-class" object
			maList <- new( "MAList", maList )

			fit <- lmFit( maList, designMatrix )
	
		} # if(!is.null(twoColour))
		
		# Otherwise this is single colour Affymetrix array, make ExpressionSet and apply linear model fit.
		else {
		
			# Subset the required columns of the expressions matrix.
			exprsForContrast <- normExprs[ , c( refAssays, testAssays ) ]
			
			# Make an ExpressionSet for this contrast.
			esetForContrast <- makeEset( exprsForContrast, refAssays, testAssays )

			# Make design matrix.
			designMatrix <- model.matrix( as.formula( "~0 + groups" ) , data = esetForContrast )
			

			# Fit linear model.
			fit <- lmFit( esetForContrast, design=designMatrix )
			
			# simple contrast matrix.
			contrast <- as.matrix( c( -1, 1 ) )

			# Apply contrast.
			fit <- contrasts.fit( fit, contrasts=contrast )
		}

		# moderated t-stats and associated p-values, B-stats (log-odds).
		fitEbayes <- eBayes( fit )

		# Independent filtering.
		# Make a data frame containing the row variances of the normalized data
		# matrix, and the unadjusted p-values.
		filterData <- data.frame( rowvar = rowVars( normExprs ), test = fitEbayes$p.value )
		
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

		# Add the column of adjusted p-values at this index to fitEbayes.
		fitEbayes$adjPvals <- filteredAdjustedPvalues[ , maxRejectionsIndex ]

		# Relevant results to data frames.
		# For 2-colour: 		
		if( !is.null( aValuesFile ) ) {
			
			# Have to specify which column of fitEbayes$t and
			# fitEbayes$coefficients we want. This is always the first column.

			# results for heatmap matrix display:
			contrastResults <- data.frame( 
				designElements = fitEbayes$genes, 
				adjPval = fitEbayes$adjPvals, 
				t = fitEbayes$t[,1], 
				logFC = fitEbayes$coefficients[,1] 
			)
			
			# results to be used for MvA plot creation (as above but with average intensities and without t-stats:
			plotData <- data.frame(
				designElements = fitEbayes$genes, 
				adjPval = fitEbayes$adjPvals, 
				logFC = fitEbayes$coefficients[,1], 
				avgExpr = fitEbayes$Amean 
			)
		}
		# For Affymetrix 1-colour:
		else {
			# results for heatmap matrix display:
			contrastResults <- data.frame(
				designElements = featureNames(esetForContrast), 
				adjPval = fitEbayes$adjPvals, 
				t = fitEbayes$t, 
				logFC = fitEbayes$coefficients
			)

			# results to be used for MvA plot creation (as above but with average intensities and without t-stats):
			plotData <- data.frame(
				designElements = featureNames(esetForContrast), 
				adjPval = fitEbayes$adjPvals, 
				logFC = fitEbayes$coefficients, 
				avgExpr = fitEbayes$Amean
			)
		}

		# write results.
		write.table(contrastResults, file=resFile, row.names=FALSE, quote=FALSE, sep="\t")
		write.table(plotData, file=plotDataFile, row.names=FALSE, quote=FALSE, sep="\t")
		
	})
}


addTechRepAverages <<- function( expressionValues, joinedTechRepAssays, twoColour ) {

	# If this is two-colour data, get the dye names: go through the joined tech
	# rep assay names and check they all share the same dye (either Cy3 or Cy5)
	# -- otherwise we can't continue. Return the dye name for each group.
	if( !is.null( twoColour ) ) {
		techRepDyes <- getDyeNames( joinedTechRepAssays )
	}

	# Get the ID (e.g. "T1") from the separator between the assay names
	# of each technical replicate group.
	techRepIDs <- gsub( "^.*<(T\\d+)_SEP>.*$", "\\1", joinedTechRepAssays )
	
	# Make technical replicate group names (hopefully no real assay names will be the same as these).
	techRepGroups <- paste( "technical_replicate_group_", techRepIDs, sep = "" )
	
	# Go through the joined assay names...
	techRepAverages <- data.frame( sapply( joinedTechRepAssays, function( x ) {
		
		# Split on the tech rep group separator, and make names R-safe.
		techRepAssayNames <- make.names( unlist( strsplit( x, "<T\\d+_SEP>" ) ) )

		# If this is a two-colour design, remove dye names from assay names.
		if( !is.null( twoColour ) ) {
			
			techRepAssayNames <- sub( ".Cy\\d$", "", techRepAssayNames )
		}
		
		# Get the normalized expression levels for these assays.
		techRepCols <- expressionValues[ , techRepAssayNames, drop=FALSE ]
		
		# Take the mean of each row.
		techRepAverages <- data.frame( apply( techRepCols, 1, function( x ) { mean( x ) } ) )
	} ) )
	
	# Name the columns of the data frame that we made using the technical replicate group names.
	colnames( techRepAverages ) <- techRepGroups

	# Now get all the technical replicate assay names for this group of assays.
	allTechRepAssayNames <- make.names( unlist( strsplit( joinedTechRepAssays, "T\\d+_SEP>" ) ) )

	# Remove these columns from expressionValues.
	normExpr_noTechReps <- expressionValues[ , !( colnames( expressionValues ) %in% allTechRepAssayNames ), drop = FALSE ]

	# Create a new data frame joining the leftover (non-tech-rep) columns and
	# the averaged (tech rep) columns.
	new_expressionValues <- cbind( normExpr_noTechReps, techRepAverages )

	# Return the new data frame.
	return( new_expressionValues )
}


replaceTechRepAssayNames <<- function( allAssayNames, twoColour ) {
		
	# Get the technical replicate assay names (joined together by e.g. "<T1_SEP>").
	joinedTechRepAssays <- allAssayNames[ grep( "<T\\d+_SEP>", allAssayNames ) ]
	
	# If this is two-colour data, get the dye names: go through the joined tech
	# rep assay names and check they all share the same dye (either Cy3 or Cy5)
	# -- otherwise we can't continue. Return the dye name for each group.
	if( !is.null( twoColour ) ) {
		techRepDyes <- getDyeNames( joinedTechRepAssays )
	}
	
	# Get the technical replicate group IDs (e.g. "T1").
	techRepIDs <- gsub( "^.*<(T\\d+)_SEP>.*$", "\\1", joinedTechRepAssays )

	# Make the technical replicate group names.
	techRepGroups <- paste( "technical_replicate_group_", techRepIDs, sep = "" )

	# Add dye names back if this is two-colour data.
	if( !is.null( twoColour ) ) {
		techRepGroups <- paste( techRepGroups, techRepDyes, sep="." )
	}
	
	# Remove the existing techincal replicate assay names from the vector of
	# assay names.
	allAssayNames <- allAssayNames[ -( grep( "<T\\d+_SEP>", allAssayNames ) ) ]

	# Add the new technical replicate group names instead.
	allAssayNames <- c( allAssayNames, techRepGroups )

	return( allAssayNames )
}


getDyeNames <<- function( joinedTechRepAssays ) {

	# Go through the groups of joined assay names...
	techRepDyes <- sapply( joinedTechRepAssays, function( x ) {
		
		# Split on tech rep group separator.
		assayNames <- unlist( strsplit( x, "<T\\d+_SEP>" ) )

		# Get a vector of all the dyes.
		dyes <- gsub( "^.*\\.(Cy\\d)$", "\\1", assayNames )

		# Check that there is exactly one dye name for this group, die if not.
		if( length( unique( dyes ) ) != 1 ) { 
			stop( "Error: did not find exactly one dye name for a two-colour technical replicate group" )
		}
		# If there is one dye name, assign to variable.
		else {
			dye <- unique( dyes )
		}
	})
	# Remove the names sapply added as we don't need them.
	names( techRepDyes ) <- NULL

	return( techRepDyes )
}


# makeEset()
# 	- Create a Biobase ExpressionSet out of the assay, contrast and expression data.
#
makeEset <<- function( exprsForContrast, refAssays, testAssays ) {

	# Expression data
	exprsForContrast <- as.matrix( exprsForContrast )
	expressionData <- assayDataNew( storage.mode = "lockedEnvironment", exprs = exprsForContrast )
	featureNames( expressionData ) <- rownames( exprsForContrast )
	sampleNames( expressionData ) <- colnames( exprsForContrast )

	# Factor data
	groups <- factor( c( rep( "ref", length( refAssays ) ), rep( "test", length( testAssays ) ) ) )
	groupsDF <- data.frame( groups = groups )
	rownames( groupsDF ) <- c( refAssays, testAssays )
	groupsDF <- new( "AnnotatedDataFrame", data = groupsDF )

	# Feature data (just design element IDs, dies without this)
	featureData <- new( "AnnotatedDataFrame", data = data.frame( designElements = rownames( exprsForContrast ) ) )
	featureNames( featureData ) <- rownames( exprsForContrast )
	
	# Return an ExpressionSet with all this data
	return( new(
		"ExpressionSet", 
		assayData = expressionData, 
		phenoData = groupsDF, 
		featureData = featureData
	) )
}


# Run with arguments if there are any, otherwise don't do anything.
args <- commandArgs( TRUE )
if( length( args ) > 0) {
	do.call( diffAtlas_DE_limma, as.list( args ) )
}
