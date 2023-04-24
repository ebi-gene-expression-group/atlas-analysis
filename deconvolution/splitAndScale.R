#!/usr/bin/env Rscript

## Split FPKMs into runs by organism part and
## Matrix scaling script for T  matrices (@zgr2788)
##

suppressMessages(library(Matrix))
suppressMessages(library(future))
suppressMessages(library(dplyr))


split_counts_per_tissue <- function(bulk_counts, design_file){
    bulk_counts_splits = list()
    # get all tissues from sdrf file that are present in experiment
    organism_parts = as.character(unique(design_file$Characteristics.organism.part.))
  
    if (length(organism_parts)== 0){   
        organism_parts = as.character(unique(design_file$Characteristics..organism.part.))
    }
    #do a check on the sdrf file   
    if (length(organism_parts)== 0){   
        stop("sdrf file does not include collum Characteristics[organism.part] or Characteristics [organism.part]. Make sure organism part names are
            stored in this collumn")
    }

    tissue_dict = c()
    for (organism_part in organism_parts){
        # get runs that belong to tissue
        tissue_runs = design_file[design_file$Characteristics.organism.part. == organism_part, ]
        tissue_runs = unique(rownames(tissue_runs))
   
        # extract column from count file for each tissue
        bulk_counts_split = bulk_counts[,colnames(bulk_counts) %in% tissue_runs]
        bulk_counts_splits <- append(bulk_counts_splits, list(bulk_counts_split))
    
    }
    names(bulk_counts_splits) = organism_parts
    return(bulk_counts_splits)   
}

args <- commandArgs(trailingOnly = TRUE)
fpkms = args[1]
sdrf_filename = args[2]
exp_name = args[3]


if(length(args) != 3) {
   stop("Not correct number of arguments. Please supply three arguments")
}

fpkms = read.csv(fpkms, sep = "\t", check.names = FALSE, row.names = 1)

# remove duplicated columns in fpkms 
fpkms = fpkms %>% 
  select(unique(colnames(.)))

sdrf = read.csv(sdrf_filename, sep = "\t", row.names = NULL)

# Rename row names to values in "Comment.ENA_RUN." or "Comment..ENA_RUN." column
if ("Comment.ENA_RUN." %in% colnames(sdrf)) {
  sdrf = sdrf[!duplicated(sdrf[ , "Comment.ENA_RUN."]),]
  rownames(sdrf) <- sdrf$Comment.ENA_RUN.
} else if ("Comment..ENA_RUN." %in% colnames(sdrf)) {
  sdrf = sdrf[!duplicated(sdrf[ , "Comment..ENA_RUN."]),]
  rownames(sdrf) <- sdrf$Comment..ENA_RUN.
} else {
  stop("Neither 'Comment.ENA_RUN.' nor 'Comment..ENA_RUN.' column found in file.")
}



# split fpkms into runs per organism part
tissue_splits = split_counts_per_tissue(fpkms, sdrf)

# set methods to scale sdrf counts
method <- 'TMM'
plan('multisession', workers = as.numeric('32')) #Paralellism
for (i in seq_along(tissue_splits)){

    # T <- readRDS(split)
    #filename_T <- sub("Tissue_splits/", "Tissue_splits/Normalized_tables/", filename_T)

    #Preprocess matrices to avoid errorshaving multiple snakefiles in one directory
    tissue = names(tissue_splits)[i]
    tissue_split = tissue_splits[[tissue]]

    message(paste0("Scaling T with method: ", method))

    #Only get rows where there is at least one cell with a count
    tissue_split <- tissue_split[rowSums(T) != 0,]

    #Same for cols
    tissue_split <- tissue_split[,colSums(tissue_split) != 0]

    #Only keep rows with different counts after transformations
    tissue_split <- tissue_split[!apply(tissue_split, 1, function(x) var(x) == 0),]

    #Switch modes
    switch(method,
    
    "TMM" = {
        suppressMessages(library(edgeR))

        tissue_split <- edgeR::DGEList(counts = tissue_split, group = colnames(tissue_split))
        tissue_split <- edgeR::calcNormFactors(tissue_split, method = "TMM")
        tissue_split <- edgeR::cpm(tissue_split)
    }
    )
    #Only get rows where there is at least one cell with a count
    tissue_split <- tissue_split[rowSums(tissue_split) != 0,]

    #Same for cols
    tissue_split <- tissue_split[,colSums(tissue_split) != 0]

    #Only keep rows with different counts after transformations
    tissue_split <- tissue_split[!apply(tissue_split, 1, function(x) var(x) == 0),]

    # Save transformed tables
    saveRDS(tissue_split, paste0('Tissue_splits/', exp_name, '/', exp_name, '-', gsub(' ', '_', tissue), '-fpkms_scaled.rds'))
}
