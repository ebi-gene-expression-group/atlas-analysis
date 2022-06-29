#!/usr/bin/env Rscript

# diffAtlas_DE_limma.R
# Microarray differential expression statistics computation for Expression Atlas.

suppressMessages( library( limma ) )
suppressMessages( library( Biobase ) )
suppressMessages( library( genefilter ) )
suppressMessages( library( reshape2 ) )

suppressMessages( library( ExpressionAtlasInternal ) )

initial.options <- commandArgs(trailingOnly = FALSE)
script.dir <- dirname(sub("--file=", "", initial.options[grep("--file", initial.options)]))
source(file.path(script.dir, 'diffAtlas_utils.R'))


# diffAtlas_DE_limma()
# - Differential expression analysis (2-group comparison) using limma.
# Arguments:
#     expAcc <- ArrayExpress accession of experiment.
 #     atlasProcessingDirectory <- path to Atlas processing directory.
diffAtlas_DE_limma <- function( expAcc, atlasProcessingDirectory ) {

   e <- try({

     # Make the config filename.
     xmlConfigFilename <- paste( expAcc, "-configuration.xml", sep = "" )
     xmlConfigFilename <- file.path( atlasProcessingDirectory, xmlConfigFilename )

     # First we need to parse the config file.
     cat( paste( "Reading XML config from", xmlConfigFilename, "..." ) )
     experimentConfig <- parseAtlasConfig( xmlConfigFilename )
     cat( "Successfully read XML config.\n" )

     # Get the experiment type.
     experimentType <- experimentConfig$experimentType

     cat( paste( "Experiment type is", experimentType, "\n" ) )

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

       cat( "Running one-colour analysis...\n" )

       run_one_colour_analysis( expAcc, allAnalytics, atlasProcessingDirectory )

     } else if( grepl( "2colour", experimentType ) ) {

       cat( "Running two-colour analysis...\n" )

       run_two_colour_analysis( expAcc, allAnalytics, atlasProcessingDirectory )

     } else {
       cat( paste(
         "Experiment type",
         experimentType,
         "not recognised. Cannot continue."
       ) )
     }

   } )

   # Die if we got an error.
   if( class( e ) == "try-error" ) {
     stop( e )
   }


}


# run_one_colour_analysis
#     - For a one-colour experiment, given an experiment accession, the list of
 #     analytics objects, and the Atlas processing directory, run the differential
 #     expression analysis defined in the contrasts contained in the analytics
 #     objects.
run_one_colour_analysis <- function( expAcc, allAnalytics, atlasProcessingDirectory ) {


   cat( paste( length( allAnalytics ), "array designs found.\n" ) )

   # Go through the analytics and do the analysis...
   invisible( sapply( allAnalytics, function( analytics ) {

     cat( paste( "Calculating differential expression statistics for array design", platform( analytics ), "...\n" ) )

     # Read experiment, contrasts and expression. Subset expression to match
     # the derived experiment.

     experiment <- exp_metadata_from_assay_groups(analytics)
     contrastsTable <- make_exp_contrast_table(analytics)

     normalizedData <- read_exp_data_table(expAcc, atlasProcessingDirectory, analytics, experiment, type = 'normalized-expressions')

     # Now we can remove duplicates (techreps) from exp
     experiment <- experiment[! duplicated(experiment$BioRepName), -1]

     cat( paste( "Found", nrow(contrastsTable), "contrasts and", nrow(experiment), "assay groups for this array design.\n" ) )

     # We'll model all the contrasts (e.g. with/without batch effects) for each formula together 

     sapply(split(contrastsTable, contrastsTable$formula), function(fc){

       # Subset to the assays relevant for this formula / group of contrasts

       normalizedDataForFormula <- normalizedData[ , colnames(normalizedData) %in% experiment$BioRepName[experiment$assay_group_id %in% c(fc$reference_assay_group_id, fc$test_assay_group_id)] ]

       cat(paste0('Processing contrasts for the "', fc$formula[1], '" formula\n'))
       designMatrix <- model.matrix(as.formula(fc$formula[1]), data=experiment)

       if( !is.fullrank( designMatrix ) ) {
         cat( "WARN  - Design matrix is not full rank, reverting to simple design matrix." )
         formulaString <- "~ groups"
         designMatrix <- model.matrix( as.formula( formulaString ), data = experiment )
       }

       # Do the first fit

       colnames(designMatrix) <- sub('assay_group_id', 'assay_group_id.', colnames(designMatrix))
       fit <- lmFit(normalizedDataForFormula, designMatrix)

       contrastNames <- paste(paste('assay_group_id', make.names(fc$test_assay_group_id), sep="."), paste('assay_group_id', make.names(fc$reference_assay_group_id), sep="."), sep="-")
       contrast.matrix <- makeContrasts(contrasts=contrastNames, levels=designMatrix)

       # Fit all the contrasts

       fit2 <- contrasts.fit(fit, contrast.matrix)
       fit2 <- eBayes(fit2)

       # This is Atlas-specific stuff:

       # Adjust the p-values and perform independent filtering, add to the fit object.

       cat( "Performing independent filtering and adjusting p-values...\n" )

       fit2$adjPvals <- do.call(cbind, lapply(1:nrow(contrastsTable), function(x){

         cat( paste( "Processing contrast", contrastsTable$contrast_name[x], "\n"))

         # Adjust the p-values and perform independent filtering
         filter_and_adjust_pvalues( rowVars( normalizedDataForFormula ), fit2$p.value[ , x ] )
       }))

       cat( "Filtering and adjustment successful.\n" )

       write_limma_de_results(expAcc, contrastsTable, fit2)

       cat( paste( "Successully completed differential expression analysis for all contrasts\n" ) )

     })

   } ) )


}


# run_two_colour_analysis
#     - For a two-colour experiment, given an experiment accession, the list of
 #     analytics objects, and the Atlas processing directory, run the differential
 #     expression analysis defined in the contrasts contained in the analytics
 #     objects. For each contrast only the samples pertinent to that contrast are
 #     modelled.
run_two_colour_analysis <- function( expAcc, allAnalytics, atlasProcessingDirectory ) {

   cat( paste( length( allAnalytics ), "array designs found.\n" ) )

   # Go through the analytics and do the analysis...
   invisible( sapply( allAnalytics, function( analytics ) {

     cat( paste( "Calculating differential expression statistics for array design", platform( analytics ), "...\n" ) )

     # Read experiment, contrasts and expression. Subset expression to match
     # the derived experiment.

     experiment <- exp_metadata_from_assay_groups(analytics, twocolor = TRUE)
     contrastsTable <- make_exp_contrast_table(analytics)

     cat( paste( "Found", nrow(contrastsTable), "contrasts and", nrow(experiment), "assay groups for this array design.\n" ) )

     # Read log fold changes and average intensities

     logFoldChanges <- read_exp_data_table(expAcc, atlasProcessingDirectory, analytics, experiment, type = 'log-fold-changes')
     averageIntensities <- read_exp_data_table(expAcc, atlasProcessingDirectory, analytics, experiment, type = 'average-intensities')

     # We'll model all the contrasts (e.g. with/without batch effects) for each formula together 

     apply(contrastsTable, 1, function(cont){

       # Die if the test and reference assay names without dye names are not identical.

       if (any( !sort(experiment$AssayNameNoCy[experiment$assay_group_id == cont['reference_assay_group_id']]) == sort( experiment$AssayNameNoCy[experiment$assay_group_id == cont['test_assay_group_id'] ]))){
         stop(  "Differing assay names found in test and reference assay groups after removing \".Cy3\" and \".Cy5\". Please verify this experiment has a two-colour design." )
       }

       # Select data specific to the formula

       contrastBioreps <- experiment$BioRepNameNoCy[experiment$assay_group_id %in% c(cont['reference_assay_group_id'], cont['test_assay_group_id'])] 

       # Do the actual modelling and analysis

       cat( "Creating targets data frame...\n" )

       targetsDF <- reshape2::dcast(subset(experiment, BioRepNameNoCy %in% contrastBioreps), BioRepNameNoCy ~ dye, value.var = 'assay_group_id')
       colnames(targetsDF)[1] <- 'BioRepName'
       rownames(targetsDF) <- targetsDF$BioRepName

       cat( "Targets data frame created successfully.\n" )

       cat( "Creating design matrix...\n" )

       designMatrix <- modelMatrix( targetsDF, ref = cont['reference_assay_group_id'] )

       cat( "Design matrix created successfully.\n" )

       cat( "Re-ordering data columns for this contrast...\n" )

       logFCsForContrast <- logFoldChanges[ , rownames( targetsDF ) ]
       avgIntsForContrast <- averageIntensities[ , rownames( targetsDF ) ]

       cat( "Data columns re-ordered successfully.\n" )

       cat( "Creating MAList object...\n" )

       maList <- list(
         genes = rownames( logFCsForContrast ),
         M = logFCsForContrast,
         A = avgIntsForContrast
       )
       maList <- new( "MAList", maList )

       cat( "MAList created successfully.\n" )

       cat( "Fitting linear model...\n" )

       fit <- lmFit( maList, designMatrix )

       cat( "Fit successful.\n" )

       cat( "Calculating differential expression statistics...\n" )

       fit <- eBayes( fit )

       cat( "Calculation successful.\n" )

       cat( "Performing independent filtering and adjusting p-values...\n" )

       # Adjust the p-values and perform independent filtering, add to the fit object.
       fit$adjPvals <- matrix(filter_and_adjust_pvalues( rowVars( logFCsForContrast ), fit$p.value[ , 1 ] ), ncol = 1)

       cat( "Filtering and adjustment successful.\n" )

       write_limma_de_results(expAcc, contrastsTable, fit)

     })
   } ))


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