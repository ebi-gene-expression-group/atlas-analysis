#!/usr/bin/env Rscript

## script to collect check correlation between the deconvolution resutls
## between FARDEEP, DWLS and EpiDISh

suppressMessages(library(tidyr))
suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)

accession <- args[1]
tissue <- args[2]

if(length(args) != 2) {
  stop("Not correct number of arguments. Please supply two arguments")
}

# in case whitespace in tissue name was replaced by '_' revert that
tissue_name <- gsub('_', ' ', tissue)
tissue_pattern <- gsub("\\(", "\\\\(", tissue)
tissue_pattern <- gsub('\\)', '\\\\)', tissue_pattern)

# check if all three deconvolution results are there
results_filenames <- paste0('Output/',accession , '/', list.files(paste0("Output/", accession), pattern=paste0(accession,'-', tissue_pattern)))
# check if all three deconvolution tools produced output
if (length(results_filenames) == 3) {
  # Read RDS files
  deconvolution_results <- lapply(results_filenames, readRDS)
  # make sure all results have the same row/colnames order
  for (i in 1:length(deconvolution_results)){
    deconvolution_results[[i]] <- deconvolution_results[[i]][,order(colnames(deconvolution_results[[i]]))]
    deconvolution_results[[i]] <- deconvolution_results[[i]][order(rownames(deconvolution_results[[i]])),]
  }
  # Calculate Pearson correlation for each RUN between the three methods
  # to see if the difference between methods is too big
  mean_vector <- numeric()
  for (i in 1:dim(deconvolution_results[[1]])[2]) {
    cor_mat <- cor(data.frame(deconvolution_results[[1]][, i], deconvolution_results[[2]][, i], deconvolution_results[[3]][, i]), method = "pearson")
    
    # Set diagonal to NA
    diag(cor_mat) <- NA
    
    mean_value <- mean(cor_mat, na.rm = TRUE)
    mean_vector[i] <- mean_value
  }
  
  # Check if mean correlation is low (0.6 seems to be a good threshold)
  if (mean(mean_vector) < 0.6) {
    cat(paste0('correlation_between_methods_lower_than_0.6'))
  } else {
    cat(paste0('mean_correlation:', round(mean(mean_vector), 3)))
  }
} else {
  cat('one_deconvolution_tool_failed')
}
