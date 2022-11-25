#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(data.table))

# Get commandline arguments.
args <- commandArgs(TRUE)

if (length(args) != 2) {
    stop("\nUsage:\n\ttranscripts_expr_values_check.R <TPMtransExpr_undecorated_filename> <RawtransExpr_undecorated_filename>\n\n")
}

# Get the transcripts expression filename from the arguments.
transExprTpmFileName <- args[1]
transExprRawFileName <- args[2]


# file checks
check_file_exists <- function(filename) {
    if (!file.exists(filename)) {
        stop(paste("Cannot find:", filename))
    }
}

check_file_exists(transExprTpmFileName)
check_file_exists(transExprRawFileName)

fread(input = transExprTpmFileName) -> TPMexpr
fread(input = transExprRawFileName) -> Rawexpr

# if any Inf exists in TPMs convert to NA
for (j in 1:ncol(TPMexpr)) set(TPMexpr, which(is.infinite(TPMexpr[[j]])), j, NA)
                 
# if any NA exists in TPMs
if (any(is.na(TPMexpr))) {

    # identify indexes of NAs in TPM transcript file
    idx = (which(is.na(TPMexpr), arr.ind = TRUE))

    # check if any associated idx in Raw transcript is 0.
    if (any(idx %in% which(Rawexpr == 0, arr.ind = TRUE))) {

        # replace NAs with 0
        TPMexpr[is.na(TPMexpr)] <- 0
        
        # if any non-numeric column
        if ( all(sapply(TPMexpr[,2:ncol(TPMexpr)], function(x){ is.numeric(x) })) == FALSE ){
            print( "ERROR - at least one column has non-numeric values" )
            quit( save = "no", status = 1, runLast = FALSE )
        }
        
        # if all zeros
        if ( all( sapply(TPMexpr[,2:ncol(TPMexpr)], function(x){ (all(x==0)  ) } )) ){
            print( "WARNING - all elements in the matrix are 0" ) 
        }

        # export TPM transcript file
        fwrite(TPMexpr, file = transExprTpmFileName, row.names = FALSE, quote = FALSE, sep = "\t")
    }
}

