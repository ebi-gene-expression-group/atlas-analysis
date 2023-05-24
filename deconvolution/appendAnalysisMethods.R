#!/usr/bin/env Rscript
### script to append analysis methods file with information (tools and references) about deconvolution analysis

suppressMessages(library('EpiDISH'))
suppressMessages(library('DWLS'))
suppressMessages(library('FARDEEP'))
suppressMessages(library(scOntoMatch))
suppressMessages(library(ontologyIndex))
suppressMessages(library('stringr'))


# Function to extract UBERON id from filename
extract_id_from_file = function(filename){
    
    base = basename(filename)
    uberon_id = str_extract(base, "UBERON_\\d{7}")
    uberon_id = sub("_", ":", uberon_id)
    return(uberon_id)
}

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

if(length(args) < 5) {
   stop("Not correct number of arguments. Please supply at least five arguments")
}

methods_file <- args[1]
accession <- args[2]
tissue <- args[3]
deconv_reference <- args[4]
workflow_base_dir <- args[5]
deconv_status <- args[6]

# read existing analysis methods file if it exits

if (!file.exists(methods_file)) {
  stop("Methods file not found.")
}
analysis_methods <- read.csv(methods_file, sep = "\t", header = F)

#get the versions of the packeages
epidish_version <- packageVersion("EpiDISH")
fardeep_version <- packageVersion("FARDEEP")
dwls_version <- packageVersion("DWLS")

# get tissue ontology
obo_file = paste0(workflow_base_dir, '/atlas-analysis/deconvolution/basic.obo') #download here: http://purl.obolibrary.org/obo/uberon/basic.obo
propagate_relationships = c('is_a', 'part_of', 'relationship: part_of', 'synonym')
# create UBRERON ontology object
ont <- ontologyIndex::get_OBO(obo_file, 
                              propagate_relationships =propagate_relationships, 
                              extract_tags = 'everything')

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
# append the version_info data frame as a new line to the analysis_methods data frame if it is not already in there
# Check if the row already exists in the DataFrame
# row_exists <- identical(which(apply(analysis_methods, 1, function(row) all(row == unlist(version_info)))), 1)

# Append the DataFrame with the new row if it doesn't already exist

if (! nrow(merge(version_info, analysis_methods))>0) {
  analysis_methods <- rbind(analysis_methods, version_info)
}

scxa_url = "<a href=https://www.ebi.ac.uk/gxa/sc/experiments/accession>accession</a>."

# add information about references 
if (deconv_status != "deconvolution_successful") {
    reference_info <- paste0(tissue, ": No suitable reference for deconvolution found")
  } else {
    deconv_tissue = getOntologyName(ont = ont, 
                                    onto_id = extract_id_from_file(deconv_reference))
    reference_info <- paste0(tissue, ": deconvolved with", deconv_tissue, "from", 
                             gsub('accession',  extract_accession(deconv_reference), scxa_url))
  }
  # create a new data frame with the reference information
  tissue_info <- data.frame(V1 = "",
                            V2 = reference_info)
                              
  
# append the tissue_info data frame as a new line to the analysis_methods data frame
analysis_methods <- rbind(analysis_methods, tissue_info)


# write the updated analysis_methods data frame to the input file with a .updated.tsv extension
write.table(analysis_methods, file = methods_file, sep = "\t", row.names = FALSE, col.names = FALSE)
