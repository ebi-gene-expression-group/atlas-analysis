#!/usr/bin/env Rscript

## script to collect devonvolution results from diffrent organism parts accross an atlas
## experiment and create final output. If deconvolution runs was not succesful for 
## an organism part, runs will not be included in file

suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))

# Function to extract accession from scxa reference
extract_accession <- function(string) {
  pattern <- "(?<=_)[A-Z]-[A-Z]+-[0-9]+"
  matches <- regmatches(string, regexpr(pattern, string, perl = TRUE))
  if (length(matches) > 0) {
    return(matches)
    } else {
    return(NA)
    }
  }

args <- commandArgs(trailingOnly = TRUE)

sdrf = args[1]
accession = args[2]
tissue =  args[3]
sc_reference = args[4]
output = args[5]
deconv_status = args[6]

if(length(args) != 6) {
   stop("Not correct number of arguments. Please supply six arguments")
}

sdrf = read.csv(sdrf, sep = "\t", row.names = NULL)
# convert "Comment.ENA_RUN." to "Comment..ENA_RUN."
colnames(sdrf) = gsub("[[:punct:]]+", '.', colnames(sdrf))
# remove rows with duplicated run ids in "Comment.ENA_RUN."
sdrf = sdrf[!duplicated(sdrf[ , "Comment.ENA.RUN."]),]
rownames(sdrf) = sdrf$Comment.ENA.RUN.
sdrf$Comment.ENA.RUN. = NULL

# in case whitespace in tissue name was replaced by '_', revert that
tissue_name = gsub('_', ' ', tissue)
# avoid problems if there is '(' or ')' in tissue name
tissue_pattern = gsub("\\(", "\\\\(", tissue)
tissue_pattern = gsub('\\)', '\\\\)', tissue_pattern)

# initialize output file if it does not exist yet
if (!file.exists(output)) {
  write.table(data.frame(ENA_RUN = character(),
                        CL_term = character(),
                        sd = numeric(),
                        proportion = numeric(),
                        organism_part = character(),
                        sc_reference = character()),
                        SCXA_experiment = character()),
                        output, sep = "\t", row.names = FALSE, col.names = TRUE)
}
# add consensus proportions to output file if deconvolution was succesfull
if (grepl('mean_correlation', deconv_status, fixed=T)){
   filenames <- paste0('Output/',accession , '/', list.files(paste0("Output/", accession), 
                                                             pattern=paste0(accession,'-', tissue_pattern)))
   
   files <- lapply(filenames, readRDS)
   # add proportions per run and cell type
   sum_props = Reduce('+',files)
   # get average of the three results
   prop = sum_props/length(files)

   # get standard deviation of the three measurements we have per run and cell type
   vec <- unlist(files, use.names = TRUE)
   DIM <- dim(files[[1]])
   n <- length(files)
   list.mean <- tapply(vec, rep(1:prod(DIM),times = n), mean)
   attr(list.mean, "dim") <- DIM
   list.mean <- as.data.frame(list.mean)
   list.sd <- tapply(vec, rep(1:prod(DIM),times = n), sd)
   attr(list.sd, "dim") <- DIM
   list.sd <- as.data.frame(list.sd)
   colnames(list.sd) = colnames(files[[1]])

   sd_norm = rowMeans(list.sd)/rowMeans(list.mean)
   rownames(list.sd) = rownames(prop)

   list.sd$CL_term = rownames(prop)

   # pivot table
   prop = as.data.frame(prop)
   prop$CL_term = rownames(prop)
   pivoted_mean = pivot_longer(prop, cols = !CL_term, names_to = "ENA_RUN", values_to = "proportion")
   pivoted_sd = pivot_longer(list.sd, cols = !CL_term, names_to = "ENA_RUN", values_to = "sd")

   MergedDF <- merge(pivoted_sd, pivoted_mean, by =c("ENA_RUN", 'CL_term'))
   MergedDF$organism_part = tissue
   MergedDF$sc_reference = sc_reference
   MergedDF$SCXA_experiment = extract_accession(basename(sc_reference))
   # append the output file with proportions
   write.table(MergedDF, output, sep = "\t", col.names = FALSE, append = TRUE, row.names = FALSE)
  } else {
   print(paste('Deconvolution not successful: no results for', tissue, 'added to', output))
  }
