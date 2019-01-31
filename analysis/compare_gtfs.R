#!/usr/bin/env Rscript

# This script compares consistency between two gtfs files 
# particularly looking percentage overalp between gene ids 
# existing in the gtfs.
# Default percentage overlap threshold : 80% 

# Load libraries
suppressMessages( library(rtracklayer) )
suppressMessages( library(optparse) )

## file checks
check_file_exists <- function( filename ) {
  if( !file.exists( filename ) ) {
    stop( paste(
      "Cannot find:",
      filename
    ) )
  }
}

## arguments
args <- parse_args(OptionParser(option_list= list(
  make_option(
    c("-e", "--ensembl_gtf"),
    help="Input gtf file from Ensembl"
  ),
  make_option(
    c("-i", "--isl_gtf"),
    help="Input gtf file from ISL"
  ),
  make_option(
    c("-p", "--percentage_threshold"),
    default = 80,
    type = 'numeric',
    help="percentage overlap"
  )
)))

## check if the file exists
check_file_exists(args$ensembl_gtf)
check_file_exists(args$isl_gtf)

# read gene ids from gtfs from ensembl and isl
gtf_file1 <- readGFF(args$ensembl_gtf, tags = c("gene_id"))
gtf_file2 <- readGFF(args$isl_gtf, tags = c("gene_id"))

## compare gene ids from gtfs files downloaded and gtf from ISL
## percentage overlap between gene ids
perc_overlap<-(length(intersect(gtf_file1$gene_id, gtf_file2$gene_id))/length(unique(gtf_file1$gene_id)))*100

## report the ones that has lower than threshold (default:80%)
if ( perc_overlap < args$percentage_threshold ) {
  ## Execution stops with exit status 1
  stop(paste ( "GTFs overlap not OK -", perc_overlap,"%","between", args$ensembl_gtf, "and", args$isl_gtf, "\n" ) )
} 
  
