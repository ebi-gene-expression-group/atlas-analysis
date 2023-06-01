#!/usr/bin/env Rscript

## script to collect check correlation between the deconvolution resutls
## between FARDEEP, DWLS and EpiDISh

suppressMessages(library(tidyr))
suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)

accession = args[1]
tissue =  args[2]

if(length(args) != 2) {
  stop("Not correct number of arguments. Please supply two arguments")
}

#in case whitespace in tissue name was replaced by '_' revert that
tissue_name = gsub('_', ' ', tissue)
tissue_pattern = gsub("\\(", "\\\\(", tissue)
tissue_pattern = gsub('\\)', '\\\\)', tissue_pattern)

# check if all three deconvolution results are there
filenames <- paste0('Output/',accession , '/', list.files(paste0("Output/", accession), pattern=paste0(accession,'-', tissue_pattern)))
# check if all three deconvolution tools produced output
if (length(filenames) == 3) {
  # Read RDS files
  files <- lapply(filenames, readRDS)
  
  # Calculate Pearson correlation for each RUN between the three methods
  # to see if the difference between methods is too big
  mean_vector <- numeric()
  for (i in 1:dim(files[[1]])[2]) {
    cor_mat <- cor(data.frame(files[[1]][, i], files[[2]][, i], files[[3]][, i]))
    
    # Set diagonal to NA
    diag(cor_mat) <- NA
    
    mean_value <- mean(cor_mat, na.rm = TRUE)
    mean_vector[i] <- mean_value
  }
  
  # Check mean correlation
  if (mean(mean_vector) < 0.7) {
    cat('correlation_too_low')
  } else {
    cat(paste0('mean_correlation:', round(mean(mean_vector), 3)))
  }
} else {
  cat('deconvolution_unsuccessful')
}
