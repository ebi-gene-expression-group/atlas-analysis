#!/usr/bin/env Rscript

# diffAtlas_DE_deseq2.R
# RNA-seq differential expression statistics computation for Expression Atlas.
#

suppressMessages( library( DESeq2 ) )

suppressMessages( library( ExpressionAtlasInternal ) )
suppressMessages( library( reshape2 ) )

initial.options <- commandArgs(trailingOnly = FALSE)
script.dir <- dirname(sub("--file=", "", initial.options[grep("--file", initial.options)]))
source(file.path(script.dir, 'diffAtlas_utils.R'))

# diffAtlas_DE_deseq2
# - Differential expression analysis (2-group comparison) using DESeq2.
# Arguments:
# 	expAcc <- ArrayExpress accession of experiment.
# 	atlasProcessingDirectory <- path to Atlas processing directory.

diffAtlas_DE_deseq2 <- function( expAcc, atlasProcessingDirectory ) {

   e <- try( {

     # Make the config filename.
     xmlConfigFilename <- paste( expAcc, "-configuration.xml", sep = "" )
     xmlConfigFilename <- file.path( atlasProcessingDirectory, xmlConfigFilename )

     # Parse the config.
     cat( paste( "Reading XML from", xmlConfigFilename, "...\n" ) )
     experimentConfig <- parseAtlasConfig( xmlConfigFilename )
     cat( "Successfully read XML config.\n" )

     # Get the experiment type.
     experimentType <- experimentConfig$experimentType

     cat( paste( "Experiment type is", experimentType, "\n" ) )

     # Check this is an RNA-seq experiment.
     if( !grepl( "rnaseq", experimentType ) ) {
       stop( paste(
         "Experiment type",
         experimentType,
         "does not look like an RNA-seq experiment. Cannot continue."
       ) )
     }

     # Get the list of analytics objects from the config.
     allAnalytics <- experimentConfig$allAnalytics

     # There should only be one analytics object for an RNA-seq experiment,
     # so make sure this is the case.
     if( length( allAnalytics ) > 1 ) {
       stop( "More than one analytics element found in XML. Cannot continue." )
     }

     # Get the analytics object.
     analytics <- allAnalytics[[ 1 ]]

     if( platform( analytics ) != "rnaseq" ) {
       stop( paste( 
         "Don't know what to do with analytics of type",
         platform( analytics )
       ) )
     }

     # Read experiment, contrasts and expression. Subset expression to match
     # the derived experiment.

     experiment <- exp_metadata_from_assay_groups(analytics)
     contrastsTable <- make_exp_contrast_table(analytics)

     # Read in the raw counts.
     cat( paste( "Reading raw counts for", expAcc, "...\n" ) )
     countsMatrix <- read_exp_data_table(expAcc, atlasProcessingDirectory, analytics, experiment, type = 'raw-counts')

     cat( "Successfully read raw counts.\n" )

     sapply(split(contrastsTable, contrastsTable$formula), function(fc){

       countsForFormula <- countsMatrix[ , colnames(countsMatrix) %in% experiment$BioRepName[experiment$assay_group_id %in% c(fc$reference_assay_group_id, fc$test_assay_group_id)] ]

       # Now we have all the info we need, we can create the DESeqDataSet object.
       cat( "Creating DESeqDataSet object...\n" )

       deseqDataSet <- DESeqDataSetFromMatrix(
         countData = countsForFormula,
         colData = experiment,
         design = as.formula( fc$formula[1] )
       )
       deseqDataSet <- DESeq( deseqDataSet )

       cat( "Differential expression analysis successful.\n" )

       cat( "Performing independent filtering and creating results table...\n" )

       formulaResults <- apply(fc, 1, function(cont){
         res <- results( deseqDataSet, contrast = c( "assay_group_id", cont['test_assay_group_id'], cont['reference_assay_group_id'] ) )

         # Make sure that all the columns of res are numeric. If not maybe
         # something went wrong; non-numeric values cause problems later.
         if( !all( sapply( res, is.numeric ) ) ) {
           stop( "Non-numeric values found in DESeq results. Cannot continue." )
         }

         res
       })

       write_deseq2_de_results(expAcc, contrastsTable, formulaResults)

       cat( "Independent filtering and results table creation successful.\n" )
     })

   } ) # try

   # Die if we got an error.

   if( class( e ) == "try-error" ) {
     stop( e )
   }
 } 

###################
# RUN MAIN FUNCTION
###################

# Run with arguments if there are any, otherwise don't do anything. Having this
# lets us source this file in an R session and load all the functions so we can
# run them separately if desired.
args <- commandArgs( TRUE )
if( length( args ) > 0) {
	do.call( diffAtlas_DE_deseq2, as.list( args ) )
}