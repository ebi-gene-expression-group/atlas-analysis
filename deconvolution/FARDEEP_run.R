#!/usr/bin/env Rscript
## script to run deconvolution with FARDEEP
##
## @zgr2788

# load FARDEEP
if (!requireNamespace("FARDEEP", quietly = TRUE)) {
  install.packages("FARDEEP", repos = "http://cran.us.r-project.org")
}

suppressMessages(library("FARDEEP"))

# get args and load files
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 3) {
   stop("Not correct number of arguments. Please supply three arguments")
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
sc_reference <- sc_reference[common, ]
fpkms  <- fpkms[common, ]
# get results and reorder the matrices for correspondence
result <- fardeep(sc_reference, fpkms, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)

p_val = result$pval
result <- t(result$abs.beta)
result <- apply(result,2,function(x) x/sum(x)) #explicit STO constraint

# save and exit
saveRDS(result, file=filename_output)
