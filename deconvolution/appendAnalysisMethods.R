#!/usr/bin/env Rscript
### script to append analysis methods file with information about deconvolution analysis

suppressMessages(library('EpiDISH'))
suppressMessages(library('DWLS'))
suppressMessages(library('FARDEEP'))

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 2) {
   stop("Not correct number of arguments. Please supply two arguments")
}

methods_file = args[1]
accession = args[2]
# read existing analysis methods file if it exits
# Load Seurat object
if (!file.exists(methods_file)) {
  stop("Methods file not found.")
}
analysis_methods <- read.csv(methods_file, sep = "\t", header = F)

#get the versions of the packeages
epidish_version <- packageVersion("EpiDISH")
fardeep_version <- packageVersion("FARDEEP")
dwls_version <- packageVersion("DWLS")

# create a new data frame with the package version information
version_info <- data.frame(
    V1 = "Deconvolution",
    V2 = paste0('Cell type proportion predictions of FPKMs per organism', 
                'part with organism part specific scRNA-seq as reference.',
                'For reference creation see',
                '<a href=https://github.com/ebi-gene-expression-group/atlas-deconvolution-references>atlas-deconvolution-references</a>.',
                'Cell type predictions displayed are averages of results from EpiDISH (version: ', as.character(epidish_version),'),',
                'FARDEEP (version: ',  as.character(fardeep_version), ') and DWLS (version: ',  as.character(dwls_version), ')')
)
# append the version_info data frame as a new line to the analysis_methods data frame
analysis_methods <- rbind(analysis_methods, version_info)

# write the updated analysis_methods data frame to the input file with a .updated.tsv extension
write.table(analysis_methods, file = paste0(accession, '-analysis-methods.tsv') , sep = "\t", row.names = FALSE, col.names = FALSE)
