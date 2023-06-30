#!/usr/bin/env Rscript

## script to find suitable reference for deconvolution 
## and return UBERON id of reference:
## 1) do zooma mapping to get UBERON id for organism part of bulk
## 2) list existing deconvolution references and try to find matching UBERON id or 
##    partent term, if no reference is found return 'noref'

# Load required packages
suppressMessages(library(devtools))
suppressMessages(library(reticulate))
suppressMessages(if (!require("scOntoMatch")) devtools::install_github("YY-SONG0718/scOntoMatch"))
suppressMessages(library(scOntoMatch))
suppressMessages(library(ontologyIndex))
suppressMessages(library(stringr))
suppressMessages(library(httr))
suppressMessages(library(jsonlite))

# Function to extract UBERON id from filename
extract_id_from_file = function(filename){
    base = basename(filename)
    uberon_id = str_extract(base, "UBERON_\\d{7}")
    uberon_id = sub("_", ":", uberon_id)
    return(uberon_id)
}

# Function to get zooma mappings for organism part as we
# run this before condensed sdrf file is created
get_semantic_tag <- function(tissue, ontology) {
  tissue <- gsub('_| ', '+', tissue)
  
  
  url <- "www.ebi.ac.uk/spot/zooma/v2/api/services/annotate?propertyValue="
  res <- GET(paste0(url, tissue) )
  stop_for_status(res)
  data <- fromJSON(rawToChar(res$content))
  #filter to only get organism part ontologies
  data = data[ grepl(ontology,(data$semanticTags) ),]
  #check confidence of zooma maping
  if (nrow(data) > 0){
    data = data[1, ]
    if (data$confidence != 'LOW'){
      semanticTag <- basename(data$semanticTags[[1]])
      return(semanticTag)
    } else {
      stop(paste('No high/good confidence mapping found for', tissue))
  }} else {
    paste('No mapping found for', tissue)
  }
    
}
                 
# Function to match tissues to a list of available tissues, if mapping is
# succesfull, mapped tissue is returned
find_reference = function(uberon_id_for_deconv, sig_dir, ont){
  
  # list fiel names in signature dir to know which tissues we have references for
  sig_files = list.files(sig_dir, '*C0_scaled.rds')
  # extract uberon ids from all filenames
  reference_uberon_ids = c()
  for (file in sig_files){
    reference_uberon_ids = append(reference_uberon_ids,  extract_id_from_file(file))
  }
  # get mapping for tissues, sort to get 'finer' terms first
  for (id in reference_uberon_ids[order(reference_uberon_ids, decreasing = T)]){
    mapping = suppressMessages(getOntoMapping(ont = ont, uberon_id_for_deconv, id))
    if (length(mapping) != 0){
     break
    }
  }
  # split into labels and ids
  tissue_to_deconvolve = names(mapping)
  mapped_to = unname(mapping)
  # if tissue could be mapped to a reference tissue return the first tissue in the list
  # not sure how else to deal with tissues that have multile possible references
  if (length(mapping) > 0 && sum(mapped_to %in% reference_uberon_ids) > 0){
    # check if term mapped to something and make sure mapped_to term is available
    reference = unname(mapping)
    reference = sub(':', '_',reference)
    return(reference[1])
  }
  else {
    return('noref')
  }
}
                 
# Run functions 
args <- commandArgs(trailingOnly = TRUE)
tissue = args[1]
deconvolution_reference_dir = args[2]
workflow_base_dir = args[3]

if(length(args) != 3) {
   stop("Not correct number of arguments. Please supply two arguments")
}

# read in UBERON ontologys
obo_file = paste0(workflow_base_dir, '/atlas-analysis/deconvolution/basic.obo') #download here: http://purl.obolibrary.org/obo/uberon/basic.obo
propagate_relationships = c('is_a', 'part_of', 'relationship: part_of', 'synonym')
# create UBRERON ontology object
ont <- ontologyIndex::get_OBO(obo_file, propagate_relationships =propagate_relationships, 
                                    extract_tags = 'everything')
# get uberon id from tissue name
uberon_id = get_semantic_tag(tissue, 'UBERON')

# check if zooma mapping was succesfull
if (length(uberon_id) < 1){
    stop(paste('no uberon id found for' ,tissue))
}

# replace  '_' in UBERON ids
uberon_id  = gsub("_", ":", uberon_id)

# check if we have a suitable sc reference to deconvolve this organism part
deconvolution_reference = find_reference(uberon_id, deconvolution_reference_dir, ont)

# print reference to sdout
cat(deconvolution_reference)
