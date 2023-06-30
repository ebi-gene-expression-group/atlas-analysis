#!/usr/bin/env Rscript

## script to run deconvolution with EpiDISH
##
## written by @zgr2788

suppressMessages(library('EpiDISH'))

# get args and load files
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Incorrect number of arguments. Please supply three arguments.")
}

filename_fpkms <- args[1]
filename_sc_reference <- args[2]
filename_output <- args[3]

fpkms <- readRDS(filename_fpkms)
sc_reference <- readRDS(filename_sc_reference)

# get intersection of genes
common <- intersect(rownames(sc_reference), rownames(fpkms))
# check if common contains values
if (length(common) == 0) {
  stop("No common genes between sc reference and fpkms found. Check that both have ENSEMBLE ids as rownames")
}
sc_reference <- sc_reference[common,]
fpkms  <- fpkms[common,]
sc_reference <- as.matrix(sc_reference)  # explicit typecasting needed

# get res and reorder the matrices for correspondence
results <- t(EpiDISH::epidish(beta.m = fpkms, ref.m = sc_reference, method = "RPC")$estF)
results <- apply(results,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
results <- apply(results, 2,function(x) x/sum(x)) #explicit STO constraint

# save and exit
saveRDS(results, file=filename_output)
