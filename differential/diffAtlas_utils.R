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
        do.call(rbind, lapply(batch_effects(cont), function(x){
          do.call(rbind, lapply(names(x@batches), function(y) data.frame(id = x@batches[[y]], variable = x@effect_name, value = y  )))
        }   ))
      }
    })))
    if (!(is.null(meta_from_contrasts))){
      assaydata <- merge(assaydata, dcast(meta_from_contrasts, id ~ variable), by.x = 'AssayName', by.y = 'id')
    }
  }
  
  rownames(assaydata) <- assaydata$AssayName

  # For two-color add separate cleaned AssayName and BioRepName columns
  
  if (twocolor){
    assaydata$AssayNameNoCy <- sub( ".Cy\\d$", "", assaydata$AssayName)
    assaydata$BioRepNameNoCy <- sub( ".Cy\\d$", "", assaydata$BioRepName)
    assaydata$dye <- sub('.*(Cy\\d+)', '\\1', assaydata$AssayName)
  }
  
  # Make sure colnames are R-safe

  colnames(assaydata) <- make.names(colnames(assaydata))

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
              paste(make.names(unlist(lapply(slot_val, function(z){ z@effect_name }))), collapse=' + ')
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

read_exp_data_table <- function(expAcc, atlasProcessingDirectory, analytics, experiment, type = NULL, dataFilename = NULL){
  
  # Get the platform (array design).
  arrayDesignPart <- ''
  if (type != 'raw-counts'){
    arrayDesignPart <- paste0('_', platform( analytics ))
  }
  
  if (is.null(dataFilename)){
    # Create the data file name from the 
    dataFilename <- paste( 
      expAcc, 
      arrayDesignPart,
      "-",
      type, 
      '.tsv.undecorated', 
      sep = "" 
    )
  
    dataFilename <- file.path( atlasProcessingDirectory, dataFilename )
  }
  
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
  
  if (type %in% c('normalized-expressions', 'raw-counts')){
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

write_limma_de_results <- function(expAcc, contrastsTable, fit2){
  
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
    resFile <- file.path( Sys.getenv( "TMPDIR" ), "tmp", paste( expAcc, contrastsTable$contrast_id[contrast_number], "analytics", "tsv", sep = "." ) )
    cat( paste( "Writing differential expression results to", resFile, "...\n" ) )
    write.table( contrastResults, file=resFile, row.names=FALSE, quote=FALSE, sep="\t" )
    cat( paste( "Results written successfully for", contrastsTable$contrast_id[contrast_number],".\n" ) )
    
    plotDataFile <- file.path( Sys.getenv( "TMPDIR" ), "tmp", paste( expAcc, contrastsTable$contrast_id[contrast_number], "plotdata", "tsv", sep = "." ) )
    cat( paste( "Writing data for MvA plot to", plotDataFile, "...\n" ) )
    write.table( plotData, file=plotDataFile, row.names=FALSE, quote=FALSE, sep="\t" )
    cat( paste("Plot data written successfully for", contrastsTable$contrast_id[contrast_number], "\n" ))
  }
}

write_deseq2_de_results <- function(expAcc, contrastsTable, results){
  
  for (contrast_number in 1:nrow(contrastsTable)){
    
    res <- results[[contrast_number]]
    
    cat( paste("Creating results data frames for", contrastsTable$contrast_id[contrast_number], "...\n" ))
    
    # Here need to select wanted columns for outut files.
    # Stats results:
    contrastResults <- data.frame(
      id = rownames( res ),
      log2FoldChange = res$log2FoldChange,
      padj = res$padj
    )
  
    # MvA plot data:
    plotData <- data.frame(
      geneID = rownames( res ),
      avgExpr = res$baseMean,
      logFC = res$log2FoldChange,
      adjPval = res$padj
    )
    
    cat( "Results data frames created successfully.\n" )
    
    # Write the files.
    resFile <- file.path( Sys.getenv( "TMPDIR" ), "tmp", paste( expAcc, contrastsTable$contrast_id[contrast_number], "analytics", "tsv", sep = "." ) )
    cat( paste( "Writing differential expression results to", resFile, "...\n" ) )
    write.table( contrastResults, file=resFile, row.names=FALSE, quote=FALSE, sep="\t" )
    cat( paste( "Results written successfully for", contrastsTable$contrast_id[contrast_number],".\n" ) )
    
    plotDataFile <- file.path( Sys.getenv( "TMPDIR" ), "tmp", paste( expAcc, contrastsTable$contrast_id[contrast_number], "plotdata", "tsv", sep = "." ) )
    cat( paste( "Writing data for MvA plot to", plotDataFile, "...\n" ) )
    write.table( plotData, file=plotDataFile, row.names=FALSE, quote=FALSE, sep="\t" )
    cat( paste("Plot data written successfully for", contrastsTable$contrast_id[contrast_number], "\n" ))
  }  
  
}