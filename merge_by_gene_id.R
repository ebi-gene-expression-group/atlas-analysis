#!/usr/bin/env Rscript

# usage: merge_by_gene_id.R outputfile fileToMerge1 ... fileToMergeN
args <- commandArgs(TRUE)

outputfile <- args[1]

if (length(args) > 2) {
  suppressMessages(library(data.table))
  first <- fread(input = args[2],
                 sep = "\t",
                 header = TRUE, key = "Gene.ID", index = "Gene.ID")

  # iterate over the rest of the files
  for (fh in args[-(1:2)]) {
    other <- fread(input = fh,
                   sep = "\t",
                   header = TRUE, key = "Gene.ID", index = "Gene.ID")
    first <- merge(first, other, all = TRUE)
  }

  fwrite(first, file = outputfile, sep = "\t", quote = FALSE, na="NA")
}
