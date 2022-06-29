#!/usr/bin/env Rscript

# Script to create and save a SimpleList object containing one or more
# Bioconductor objects storing expression data from an Expression Atlas
# experiment. These object are SummarizedExperiment (RNA-seq data),
# ExpressionSet (1-colour microarray data), or MAList (2-colour microarray
# data).

suppressMessages(library(ExpressionAtlasInternal))
suppressMessages(library(optparse))

args <- parse_args(OptionParser(option_list = list(
  make_option(
    c("-s", "--source"),
    help = "Source folder"
  ),
  make_option(
    c("-a", "--accession"),
    help = "Experiment accession"
  ),
  make_option(
    c("-o", "--output"),
    help = "Output file destination"
  )
)))

print('Searching $SDRF_PATH environment variable:')
print(system('echo $SDRF_PATH', intern = TRUE))

# Create experiment summary SimpleList.
experiment_summary <- summarizeAtlasExperiment(
  args$accession,
  args$source
)

# Save the object to a file.
save(experiment_summary, file = args$output)