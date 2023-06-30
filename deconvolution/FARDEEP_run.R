#!/usr/bin/env Rscript
## script to run deconvolution with FARDEEP
##
## @zgr2788

#Load FARDEEP
if (!require("FARDEEP", quietly = TRUE))
    install.packages("FARDEEP", repos='http://cran.us.r-project.org')
suppressMessages(library("FARDEEP"))

#Get args and load files
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 3) {
   stop("Not correct number of arguments. Please supply three arguments")
}

filename_fpkms <- args[1]
filename_sc_reference <- args[2]
filename_output <- args[3]

fpkms <- readRDS(filename_fpkms)
sc_reference <- readRDS(filename_sc_reference)

#Toss out the genes tossed out in T normalization from C as well
common <- intersect(rownames(sc_reference), rownames(fpkms))
sc_reference <- sc_reference[common, ]
fpkms  <- fpkms[common, ]
#Get results and reorder the matrices for correspondence
print(length(common)) # remove later
result <- fardeep(sc_reference, fpkms, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)

p_val = result$pval
result <- t(result$abs.beta)
result <- apply(result,2,function(x) x/sum(x)) #explicit STO constraint
#if (filename_P != 'Modules/Psuedobulk/dummy_props.rds') res <- res[order(match(rownames(res), rownames(P))),]

# save and exit
saveRDS(result, file=filename_output)
