#!/usr/bin/env Rscript
## script to run deconvolution with DWLS
##
## written by @zgr2788 

# Load DWLS
suppressMessages(library(future.apply))
suppressMessages(library(DWLS))

# get args and load files
args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 6) {
   stop("Not correct number of arguments. Please supply six arguments to DWLS")
}

filename_fpkms <- args[1]
filename_sc_reference <- args[2]
filename_phenData <- args[3]
plan('multisession', workers = as.numeric(args[4])) #Paralellism
filename_output <- args[5]
scratchDir <- args[6]

fpkms <- readRDS(filename_fpkms)
sc_reference <- readRDS(filename_sc_reference)
phenData <- readRDS(filename_phenData)

# match genes in rows for both references
common <- intersect(rownames(sc_reference), rownames(fpkms))
# check if common contains values
if (length(common) == 0) {
  stop("No common genes between sc reference and fpkms found. Check that both have ENSEMBLE ids as rownames")
}
sc_reference <- sc_reference[common,]
fpkms <- fpkms[common,]
rm(common)

# preprocess
message("Started running sig...")
Signature <- buildSignatureMatrixMAST(scdata = sc_reference, 
                                      id = phenData[,"cellType"], 
                                      path = scratchDir, 
                                      diff.cutoff = 0.5, 
                                      pval.cutoff = 0.01)

# get results and reorder the matrices for correspondence
result <- future_apply(fpkms, 2, function(x) {
  b <- setNames(x, rownames(fpkms))
  tr <- trimData(Signature, b)
  RES <- t(solveDampenedWLS(tr$sig, tr$bulk))
}, future.seed = TRUE)

rownames(result) <- as.character(unique(phenData$cellType))
# convert values smaller than 0.00001 t0 0
result[result < 10^-5] <- 0 #Convergence error tolerance = 10^-5

# save and exit
saveRDS(result, file=filename_output)
print('DWLS run succesfully completed.')
