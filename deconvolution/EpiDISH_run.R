#!/usr/bin/env Rscript

## script to run deconvolution with EpiDISH
##
## @zgr2788

suppressMessages(library('EpiDISH'))

# Get args and load files
args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 3) {
   stop("Not correct number of arguments. Please supply three arguments")
}

filename_fpkms <- args[1]
filename_sc_reference <- args[2]
filename_O <- args[3]

fpkms <- readRDS(filename_fpkms)
sc_reference <- readRDS(filename_sc_reference)


# Toss out the genes tossed out in T normalization from C as well
common <- intersect(rownames(sc_reference), rownames(fpkms))
sc_reference <- sc_reference[common,]
fpkms  <- fpkms[common,]
sc_reference <- as.matrix(sc_reference) #Explicit typecasting needed

# Get res and reorder the matrices for correspondence
results <- t(EpiDISH::epidish(beta.m = fpkms, ref.m = sc_reference, method = "RPC")$estF)
results <- apply(results,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
results <- apply(results, 2,function(x) x/sum(x)) #explicit STO constraint

# Save and exit
saveRDS(results, file=filename_O)
