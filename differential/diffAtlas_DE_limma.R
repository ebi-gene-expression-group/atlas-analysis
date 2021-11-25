#!/usr/bin/env Rscript

# diffAtlas_DE_limma.R
# Microarray differential expression statistics computation for Expression Atlas.

suppressMessages( library( limma ) )
suppressMessages( library( Biobase ) )
suppressMessages( library( genefilter ) )
suppressMessages( library( reshape2 ) )

suppressMessages( library( ExpressionAtlasInternal ) )

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

# make_assays_to_bioreps_df
#     - Given a list of biological replicate objects and optional twoColour flag,
#     create a data frame mapping assay names to their biological replicate
#     names. These will either be the same as the assay name (no technical
#     replicates), or the technical replicate group ID (technical replicates).
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

# filter_and_adjust_pvalues
#     - Given row variances of the expression data, and the raw (unadjusted)
#     p-values, run independent filtering via genefilter package. Return vector
#     of BH-adjusted p-values.
filter_and_adjust_pvalues <- function( normDataRowVars, rawPvalues ) {

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
        filter = filterData$rowvar,    # use the row variances as the filter statistic.
        test = filterData$test,    # the unadjusted p-values.
        theta = theta,    # the range of filtering proportions.
        method = "BH"    # use Benjamini-Hochberg for p-value adjustment.
    )

    # filteredAdjustedPvalues is a matrix of adjusted p-values, with a
    # column for each proportion from theta.  For each column, count how
    # many rejections of the null hypothesis we will make, using the FDR
    # threshold of 0.05.
    numRej <- colSums( filteredAdjustedPvalues < 0.05, na.rm=TRUE )
    
    # Find the index of the column that had the most rejections of the null
    # hypothesis.
    maxRejectionsIndex <- which.max( numRej )

    # Return the column of adjusted p-values at this index to the fit.
    return( filteredAdjustedPvalues[ , maxRejectionsIndex ] )
}

# Just make a simple metadata table, including any batch effects  

exp_metadata_from_assay_groups <- function(analytics, twocolor = FALSE){
  
  ags <- assay_groups( analytics )
  contrasts <- atlas_contrasts( analytics )
  
  metaRows <-lapply(ags, function(x){
    df <- make_assays_to_bioreps_df(x@biological_replicates)
    
    # Add assay group info to all assays
    
    for (slot_name in slotNames(x)[slotNames(x) != 'biological_replicates']){
      slot_vals <- slot(x, slot_name)
      if (length(slot_vals) > 0){
        df[[slot_name]] <- slot(x, slot_name)
      }
    }
    df
  })
  
  # Allow for assay groups with different variables
  
  all_colnames <- unique(unlist(lapply(metaRows, function(x) names(x))))
  assaydata <- data.frame(do.call(rbind, lapply(metaRows, function(x){ x[! names(x) %in% all_colnames] <- NA; x[all_colnames]})))
  
  # Subset assay data to only the groups we'll be making contrasts from
  
  all_contrast_assaygroups <- unique(unlist(lapply(contrasts, function(x) c(x@reference_assay_group_id, x@test_assay_group_id))))
  assaydata <- assaydata[assaydata$assay_group_id %in% all_contrast_assaygroups,]
  
  # Get covariates from contrasts (would it really be so hard to store it next to the assays?)
  
  if (any(unlist(lapply(contrasts, function(x) length(batch_effects(x)) > 0)))){
    meta_from_contrasts <- unique(do.call(rbind, lapply(contrasts, function(cont){
      if (length(batch_effects(cont)) > 0){
        if (twocolor){
          stop( "Cannot handle batch effects for 2-colour microarray data." )
        }
        do.call(rbind, lapply(cont@batch_effects, function(x){
          do.call(rbind, lapply(names(x@batches), function(y) data.frame(id = x@batches[[y]], variable = x@effect_name, value = y  )))
        }   ))
      }
    })))
    if (!(is.null(meta_from_contrasts))){
      assaydata <- merge(assaydata, dcast(meta_from_contrasts, id ~ variable), by.x = 'AssayName', by.y = 'id')
    }
  }
  
  # For two-color add separate cleaned AssayName and BioRepName columns
  
  if (twocolor){
    assaydata$AssayNameNoCy <- sub( ".Cy\\d$", "", assaydata$AssayName)
    assaydata$BioRepNameNoCy <- sub( ".Cy\\d$", "", assaydata$BioRepName)
    assaydata$dye <- sub('.*(Cy\\d+)', '\\1', assaydata$AssayName)
  }
  
  assaydata
}

# Make a simple contrasts table

make_exp_contrast_table <- function(analytics){
  
  conts <- atlas_contrasts( analytics )
  
  conts <- data.frame(do.call(
    rbind, 
    lapply(conts, function(x){ 
      unlist(
        structure(
          lapply(slotNames(x), function(y){ 
            slot_val <- slot(x,y)
            if (y == 'batch_effects' && length(slot_val) > 0){
              paste(unlist(lapply(slot_val, function(z){ z@effect_name })), collapse=' + ')
            }else{
              slot_val
            }
          }), 
          names = slotNames(x))
      )
    })
  ), stringsAsFactors = FALSE)
  
  # Make the formula we need for each contrast based on any batch effects present
  
  conts$formula <- apply(conts, 1, function(x){
    formula <- '~ 0 + assay_group_id'
    if ('batch_effects' %in% names(x) && !is.na(x['batch_effects'])){
      formula <- paste(formula, x['batch_effects'], sep = ' + ')
    }
    formula
  })
  
  conts
}

# Read a normalised expression table, log fold change table or mean intensities table

read_exp_data_table <- function(expAcc, atlasProcessingDirectory, analytics, experiment, type){
  
  # Get the platform (array design).
  arrayDesign <- platform( analytics )
  
  # Create the normalized data file name.
  dataFilename <- paste( 
    expAcc, 
    "_", 
    arrayDesign, 
    "-",
    type, 
    '.tsv.undecorated', 
    sep = "" 
  )
  
 dataFilename <- file.path( atlasProcessingDirectory,dataFilename )
  
  if( !file.exists( dataFilename ) ) {
    stop( paste( "Cannot find:", dataFilename ) )
  }
  
  cat( paste( "Reading", type, "data from", dataFilename, "...\n" ) )
  
  # Read in the normalized data.
  parsedData <- read.delim( 
    dataFilename, 
    header = TRUE, 
    stringsAsFactors = FALSE,
    row.names = 1
  )
  
  cat( paste("Successfully read", type, "data.\n") )
  cat( "Matching data to experiment.\n" )
  
  if (type == 'normalized-expressions'){
    assayNames <- experiment$AssayName
    bioRepNames <- experiment$BioRepName
  }else{
    
    # For two color (retrieving fold changes and mean intensities), we have to
    # account for the the Cy3/5 for each assay. 
    
    assayToBiorep <- unique(experiment[,c('AssayNameNoCy', 'BioRepNameNoCy')])
    
    assayNames <- as.character(assayToBiorep[,'AssayNameNoCy'])
    bioRepNames <- as.character(assayToBiorep[,'BioRepNameNoCy'])
  }
  
  parsedData <- parsedData[, assayNames]
  
  # Merge technical replicates by mean
  
  if (any(duplicated(bioRepNames))){
    parsedData <- do.call(cbind, lapply(split(data.frame(t(parsedData), check.names = FALSE), bioRepNames), colMeans))
  }
  
  parsedData
}

# Write d/e results to files

write_de_results <- function(expAcc, contrastsTable, fit2){
  
  for (contrast_number in 1:nrow(contrastsTable)){
    
    cat( paste("Creating results data frames for", contrastsTable$contrast_id[contrast_number], "...\n" ))

    contrastResults <- data.frame(
      designElements = rownames(fit2$p.value), 
      adjPval = fit2$adjPvals[,contrast_number], 
      t = fit2$t[,contrast_number], 
      logFC = fit2$coefficients[ , contrast_number ]
    )
    
    plotData <- data.frame(
      designElements = rownames(fit2$p.value), 
      adjPval = fit2$adjPvals[,contrast_number], 
      logFC = fit2$coefficients[ , contrast_number ], 
      avgExpr = fit2$Amean
    )
    cat( "Results data frames created successfully.\n" )
    
    # Write the files.
    resFile <- file.path( Sys.getenv( "HOME" ), "tmp", paste( expAcc, contrastsTable$contrast_id[contrast_number], "analytics", "tsv", sep = "." ) )
    cat( paste( "Writing differential expression results to", resFile, "...\n" ) )
    write.table( contrastResults, file=resFile, row.names=FALSE, quote=FALSE, sep="\t" )
    cat( paste( "Results written successfully for", contrastsTable$contrast_id[contrast_number],".\n" ) )
    
    plotDataFile <- file.path( Sys.getenv( "HOME" ), "tmp", paste( expAcc, contrastsTable$contrast_id[contrast_number], "plotdata", "tsv", sep = "." ) )
    cat( paste( "Writing data for MvA plot to", plotDataFile, "...\n" ) )
    write.table( plotData, file=plotDataFile, row.names=FALSE, quote=FALSE, sep="\t" )
    cat( paste("Plot data written successfully for", contrastsTable$contrast_id[contrast_number], "\n" ))
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
    
    normalizedData <- read_exp_data_table(expAcc, atlasProcessingDirectory, analytics, experiment, 'normalized-expressions')
    
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
      
      write_de_results(expAcc, contrastsTable, fit2)
      
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
    
    logFoldChanges <- read_exp_data_table(expAcc, atlasProcessingDirectory, analytics, experiment, 'log-fold-changes')
    averageIntensities <- read_exp_data_table(expAcc, atlasProcessingDirectory, analytics, experiment, 'average-intensities')
    
    # We'll model all the contrasts (e.g. with/without batch effects) for each formula together 
    
    apply(contrastsTable, 1, function(cont){
      
      saveRDS(cont, file="cont.rds")
      
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
      
      write_de_results(expAcc, contrastsTable, fit)
      
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
