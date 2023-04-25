#!/usr/bin/env Rscript

### script to create file with information about which reference
### was used for deconvolution per organism part

#Read data
args <- commandArgs(trailingOnly = TRUE)
tissue <- args[1]
accession <- args[2]
deconv_reference <- args[3]

# Define file name and data to be appended
file_name <- paste0(accession, "-deconvolution_info.tsv")
new_data <- data.frame(
  task_completed = "YES",
  accession = accession,
  tissue = tissue,
  deconv_reference = basename(deconv_reference)
  #time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

# Check if file exists
if(file.exists(file_name)) {
  # Read in existing data and check for duplicates
  existing_data <- read.table(file_name, header = TRUE, sep = "\t")
  new_rows <- !duplicated(rbind(existing_data, new_data))
  
  # Append new data to existing file if it's not a duplicate
  if(sum(new_rows) > nrow(existing_data)) {
    all_data <- rbind(new_data)
    write.table(all_data, file = file_name, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
  } else {
    cat("Duplicate row detected. No data added.\n")
  }
} else {
  # Create new file with header and new data
  write.table(new_data, file = file_name, sep = "\t", row.names = FALSE, col.names = TRUE, append = FALSE)
}
