#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(data.table))

# Get commandline arguments.
args <- commandArgs( TRUE )

if( length( args ) != 2 ) {
  stop( "\nUsage:\n\ttranscripts_expr_values_check.R <TPMtransExpr_undecorated_filename> <RawtransExpr_undecorated_filename>\n\n" )
}

# Get the transcripts expression filename from the arguments.
transExprTpmFileName <- args[ 1 ]
transExprRawFileName <- args[ 2 ]


## file checks
check_file_exists <- function( filename ) {
  if( !file.exists( filename ) ) {
    stop( paste( "Cannot find:", filename ) )
  }
}

check_file_exists(transExprTpmFileName)
check_file_exists(transExprRawFileName)

fread(input=transExprTpmFileName)->TPMexpr
fread(input=transExprRawFileName)->Rawexpr

# if any NA exists in TPMs
if ( any(is.na(TPMexpr)) ) {
  
  #identify indexes of NAs in TPM transcript file
  idx=(which(is.na(TPMexpr), arr.ind=TRUE))
  
  # check if any associated idx in Raw transcript is 0.
  if (any(idx %in% which(Rawexpr==0, arr.ind=TRUE))) {
    
    # replace NAs with 0
    TPMexpr[is.na(TPMexpr)] <- 0
    
    # export TPM transcript file
    fwrite( TPMexpr, file = transExprTpmFileName, row.names=FALSE, quote=FALSE, sep="\t" )
  }
}
