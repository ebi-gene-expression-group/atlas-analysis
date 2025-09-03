
# Identification of markers genes for baseline proteomics experiments in Expression Atlas.


required_packages <- c("MGFR", "xml2", "dplyr", "data.table") 

# function to check and load a package, fail if not installed
require_or_fail <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(paste("Package", pkg, "is required but not installed."), call. = FALSE)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# check and load each
for (pkg in required_packages) {
    require_or_fail(pkg)
}
cat("All required packages are loaded successfully.\n")



 # parse arguments
args <- commandArgs(trailingOnly = FALSE)

# Get path to the currently running script
script.path <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
script.dir <- dirname(script.path)

source(file.path(script.dir, "marker_gene_utils.R"))

args <- args[-(1:5)]
print(args)
 
 if (length(args) < 4 || length(args) > 6) {
        stop("Usage: Rscript script.R <config_xml> <expression_data.undecorated.aggregated> <expression_data> <output_file> [<specificity_score_cutoff> <expression_cutoff>]")
 }
 
 # Default values
 default_specificity_cutoff <- 0.25  # between 0 and 1
 default_expression_cutoff  <- 0.5
 
 if (length(args) == 4){
     args[5] <- default_specificity_cutoff
     args[6] <- default_expression_cutoff
 }
 
 if (length(args) == 5){
     args[6] <- default_expression_cutoff
 }


arg_names <- c("Configuration XML file", "Expression data file (undecorated)", "Expression data file", "Output file", "Specificity Score Cutt-off", "Expression Cutoff")
names(args) <- arg_names

cat("Arguments received:\n")
for (name in names(args)) {
    cat(name, ":", args[[name]], "\n")
}
cat("\n")

SPECIFICITY_SCORE_CUTOFF  <- args[5]
EXPRESSION_CUTOFF  <- args[6]

metric <- 'ppb'
accession <- as.character( sub("\\.tsv\\.undecorated.aggregated$", "", args[2]) )



####################################
# read xml file
####################################
doc <- read_xml(args[1])

# find all assay_group nodes
groups <- xml_find_all(doc, ".//assay_group")

# extract information into a data frame
assay_df <- bind_rows(lapply(groups, function(group) {
    group_id <- xml_attr(group, "id")
    group_label <- xml_attr(group, "label")
    assays <- xml_find_all(group, ".//assay")
    
    data.frame(
        group = group_id,
        label = group_label,
        assay = xml_text(assays),
        tech_rep = xml_attr(assays, "technical_replicate_id"),
        stringsAsFactors = FALSE
    )
}))

print(assay_df)

# summary per group (number of assays regardless of technical replicate info)
summary_df <- assay_df %>%
    group_by(group, label) %>%
    summarise(n_assays = n(), .groups = "drop")

print(summary_df)




################################################
# Read <accession>.tsv.undecorated.aggregated
################################################

dt <- fread( args[2] )

# Ensure dt is a data.table
setDT(dt)

# dt columns names in the form: g1.WithInSampleAbundance g2.WithInSampleAbundance etc

# Keep only 'Gene ID' and all columns with 'WithInSampleAbundance' (case-insensitive)
dt_temp <- dt[, c("Gene ID", grep("WithInSampleAbundance", names(dt), ignore.case = TRUE, value = TRUE)), with = FALSE]
dt <- dt_temp

# Keep rows where Gene ID is not NA, not empty, and not purely numeric
dt <- dt[!is.na(`Gene ID`) & `Gene ID` != "" & !grepl("^\\d+$", `Gene ID`)]

# rename columns to assay name
setnames(dt, old = names(dt)[-1], new = sapply(names(dt)[-1], function(n) assay_df$assay[assay_df$group == sub("\\..*", "", n)][1]))


# remove any possible duplicated columns with same name
dt <- dt[, !duplicated(names(dt)), with = FALSE]


# Create a named vector: names are old (assay), values are new (label)
assay_to_label <- setNames(assay_df$label, assay_df$assay)

# Get intersecting assays (columns in dt that match assays in assay_df)
columns_to_rename <- intersect(names(assay_to_label), colnames(dt))

# Rename those columns
setnames(dt, columns_to_rename, assay_to_label[columns_to_rename])

# remove columns not included in the config XML
dt <- dt[, c(1, which(names(dt)[-1] %in% assay_df$label ) + 1), with = FALSE]


# Exclude "Gene ID" column
dt_numeric <- dt[, -1, with = FALSE]
gene_ids <- dt[[1]]  # save Gene ID separately - In proteomics, sometimes this is actually the gene name

# get the column names
col_names <- colnames(dt_numeric)
unique_names <- unique(col_names)

# For each unique name, average all columns with that name
averaged_list <- lapply(unique_names, function(name) {
    cols <- which(col_names == name)
    rowMeans(dt_numeric[, ..cols], na.rm = TRUE)
})

# Combine into a new data.table with Gene ID
result_dt <- data.table(`Gene ID` = gene_ids, setNames(as.data.table(averaged_list), unique_names))


# Replaces NA with 0. Skips the first column (Gene ID)
result_dt[, (2:ncol(result_dt)) := lapply(.SD, function(x) fifelse(is.na(x), 0, x)), .SDcols = 2:ncol(result_dt)]


# Assuming 'Gene ID' is the first column and the rest are numeric
# remove Genes with no expression in any column (won't be used for marker identification)
result_dt <- result_dt[rowSums(result_dt[, -1, with = FALSE] != 0) > 0]


# Remove 'Gene ID' column and convert remaining data to matrix
mat <- as.matrix(result_dt[, !"Gene ID", with = FALSE])

# Set rownames from 'Gene ID'
rownames(mat) <- result_dt[["Gene ID"]]

# Remove rows with NA row names or row names that are purely numeric
mat <- mat[!is.na(rownames(mat)) & !grepl("^\\d+$", rownames(mat)), ]



####################################
# run MGFR to get marker genes
####################################

# Function to detect marker genes normalised expression data
markers.list <- getMarkerGenes.rnaseq(mat, 
                                      class.vec = colnames(mat),
                                      samples2compare="all", 
                                      annotate=FALSE, 
                                      gene.ids.type="ensembl", 
                                      score.cutoff=1)

print( names(markers.list) )

############################################################################
# Sanity check -  that marker genes are not repeated in any list
############################################################################

# Extract gene IDs from each list and clean them
gene_lists <- lapply(markers.list, function(x) sub(" :.*", "", x))

# combine all gene IDs into a single vector with list (tissue) names
gene_with_source <- unlist(gene_lists, use.names = FALSE)
sources <- rep(names(gene_lists), lengths(gene_lists))

# create a data frame of gene ID and source tissue
df <- data.frame(gene = gene_with_source, source = sources, stringsAsFactors = FALSE)

# count in how many distinct groups each gene appears

duplicated_genes <- df %>%
    distinct(gene, source) %>%        # remove within-list duplicates
    group_by(gene) %>%
    summarise(n_sources = n()) %>%
    filter(n_sources > 1)

# show genes present in more than one list
print(duplicated_genes)

# Stop execution if duplicates are found
if (nrow(duplicated_genes) > 0) stop("Error: Duplicated genes found across gene marker lists.")


############################################################################
# Filter all marker lists by SPECIFICITY_SCORE_CUTOFF
############################################################################
markers.list <- lapply(markers.list, function(marker_vec) {
    # Keep only entries that contain ":"
    marker_vec <- marker_vec[grepl(":", marker_vec)]
    
    # Extract scores
    scores <- as.numeric(sub(".*: ", "", marker_vec))
    
    # Filter out entries with score >= cutoff
    marker_vec[!is.na(scores) & scores <= SPECIFICITY_SCORE_CUTOFF]
})

total_markers <- sum(sapply(markers.list, length))
if (total_markers == 0) {
    message <- "No markers found after applying specificity score cutoff. Consider lowering the cutoff value."
    info_file <- paste0(args[4], ".info")
    writeLines(message, info_file)
    stop(message, call. = FALSE)
}

print( markers.list )


############################################################################
# save table to data frame
############################################################################

# initialize list to store results
marker_tables <- list()

# Loop over each entry in markers.list
for (name in names(markers.list)) {
    marker_vec <- markers.list[[name]]
    
    # Skip if fewer than 1 markers
    if (length(marker_vec) < 1) next
    
    # Split into gene and score
    split_vals <- strsplit(marker_vec, " : ")
    GENE_ID <- sapply(split_vals, `[`, 1)
    specificity_score <- as.numeric(sapply(split_vals, `[`, 2))
    
    df <- data.frame(
        ACCESSION=accession ,
        GROUP= summary_df %>%  filter(label == sub("_markers$", "", name) ) %>% pull(group),
        GROUP_NAME=sub("_markers$", "", name),
        GENE_ID = GENE_ID,
        SPECIFICITY_SCORE = specificity_score,
        RANKING = seq_along(GENE_ID),
        METRIC = metric,
        NUMBER_SAMPLES= summary_df %>%  filter(label == sub("_markers$", "", name) ) %>% pull(n_assays),
        EXPRESSION = 0
    )
    
    marker_tables[[name]] <- df
}


# Combine all data frames into one - List of final expression markers
all_markers_df <- do.call(rbind, marker_tables)   

# view result
print( head(all_markers_df ) )
print( dim(all_markers_df ) )


#############################################################################################################
# update the table to include rows for all samples, for all these genes identified as marker in one sample  
#############################################################################################################

# extract required columns from summary_df
tissue_map <- summary_df %>% select(GROUP = group, GROUP_NAME = label, NUMBER_SAMPLES = n_assays)

# Get all unique GENE_IDs and existing combinations
unique_genes <- unique(all_markers_df$GENE_ID)

# Create all possible gene-tissue combinations
complete_combos <- expand.grid(GENE_ID = unique_genes, GROUP_NAME = tissue_map$GROUP_NAME, stringsAsFactors = FALSE) %>%
    left_join(tissue_map, by = "GROUP_NAME")  # add GROUP info

# Identify missing combinations
existing_keys <- all_markers_df %>% select(GENE_ID, GROUP_NAME, GROUP, NUMBER_SAMPLES)
missing_combos <- anti_join(complete_combos, existing_keys, by = c("GENE_ID", "GROUP_NAME", "GROUP", "NUMBER_SAMPLES"))

# create placeholder rows for missing combinations
filled_rows <- missing_combos %>%
    mutate(
        SPECIFICITY_SCORE = -1,
        RANKING = -1,
        METRIC =  metric,
        EXPRESSION = 0,
        ACCESSION = accession
    )

# bind with original data
final_df <- bind_rows(all_markers_df, filled_rows)

# arrange for neatness
final_df <- final_df %>%
    arrange( GENE_ID, GROUP_NAME )


##############################################################################
# add correct gene expression
# from <accession>.tsv 
##############################################################################

expression_decorated <- read.csv2(args[3], sep="\t")

print( head(expression_decorated) )


name_map <- setNames(summary_df$label, summary_df$group)

# Rename expression_decorated columns
new_names <- colnames(expression_decorated)

# Loop over group codes and replace them in column names
for (code in names(name_map)) {
    new_names <- sub(code, name_map[[code]], new_names)
}

# if necessary: remove '.WithInSampleAbundance' from the column names
new_names <- sub("\\.WithInSampleAbundance", "", new_names)

# Assign cleaned column names
colnames(expression_decorated) <- new_names



# add column with gene_name, join on gene ID
final_df <- final_df %>%
    left_join(expression_decorated %>% select(GENE_ID = Gene.ID, GENE_NAME = Gene.Name),
              by = "GENE_ID")



# Ensure EXPRESSION is numeric
final_df$EXPRESSION <- NA_real_

# Loop over rows to update EXPRESSION based on matching GENE_ID and GROUP
for (i in seq_len(nrow(final_df)) ) {
    gene_id <- final_df$GENE_ID[i]
    group <- final_df$GROUP_NAME[i]
    
    # Check if gene exists in expression_decorated and group column exists
    row_match <- which(expression_decorated$Gene.ID == gene_id)

    if (length(row_match) == 1 && group %in% colnames(expression_decorated)) {

        expr_values <- as.numeric(expression_decorated[row_match, group]) # as.numeric(unlist(strsplit(expr_string, ",")))

        final_df$EXPRESSION[i] <- expr_values
    }
}



###############################################################
# apply cutoff of gene expression and reorder ranking
###############################################################

# Identify genes where any row has RANKING != -1 AND EXPRESSION <= EXPRESSION_CUTOFF
genes_to_remove <- final_df %>%
    filter(RANKING != -1, EXPRESSION <= EXPRESSION_CUTOFF) %>%
    distinct(GENE_ID) %>%
    pull(GENE_ID)

# remove all rows for those genes
filtered_df <- final_df %>%
    filter(!GENE_ID %in% genes_to_remove)

# reassign RANKING within each GROUP_NAME for remaining valid rows
# Re-rank genes within each group by ascending specificity, then descending expression, then gene ID (if RANKING != -1)

filtered_df <- filtered_df %>%
  group_by(GROUP_NAME) %>%
  arrange(SPECIFICITY_SCORE, desc(EXPRESSION), GENE_ID, .by_group = TRUE) %>%
  mutate(
    new_rank = NA_integer_,  # initialize
    new_rank = replace(new_rank, RANKING != -1, row_number()),  # assign ranks where applicable
    RANKING = if_else(is.na(new_rank), -1L, new_rank)
  ) %>%
  select(-new_rank) %>%
  ungroup()



# get unique genes ordered by RANKING within each group using reframe()
ordered_genes <- filtered_df %>%
    filter(RANKING > 0) %>%
    group_by(GROUP_NAME) %>%
    arrange(RANKING, .by_group = TRUE) %>%
    reframe(GENE_ID = unique(GENE_ID)) %>%
    pull(GENE_ID)


# Use those ordered genes to sort the full dataframe
filtered_df_sorted <- filtered_df %>%
    filter(GENE_ID %in% ordered_genes) %>%
    mutate(GENE_ID = factor(GENE_ID, levels = ordered_genes)) %>%
    arrange(GENE_ID, GROUP_NAME)


# replace -1 with NULL for sql loading
filtered_df_sorted <- filtered_df_sorted %>%
    mutate(
        RANKING = ifelse(RANKING == -1, "NULL", as.character(RANKING)),
        SPECIFICITY_SCORE = ifelse(SPECIFICITY_SCORE == -1, "NULL", as.character(SPECIFICITY_SCORE))
    )


colnames(filtered_df_sorted) <- c("experiment_accession", "assay_id", "assay", "gene_id", "specificity_score", "marker_gene_rank", "expression_unit", "number_assays", "expression_level", "gene_name")

# remove 'assay_id' and 'number_assays' columns, for now
filtered_df_sorted <- filtered_df_sorted[, !(names(filtered_df_sorted) %in% c("number_assays"))]

# replace NA in expression_level with 0
filtered_df_sorted$expression_level[is.na(filtered_df_sorted$expression_level)] <- 0

print ( head(filtered_df_sorted) )

filtered_df_sorted_resorted <- filtered_df_sorted[c(setdiff(names(filtered_df_sorted), "assay_id"), "assay_id")]

filtered_df_sorted <- filtered_df_sorted_resorted
                     
print ( head(filtered_df_sorted_resorted) )

# temporal fix - sometimes and Gene.Name column is empty (Gene.ID is the Gene Name)
if (all(is.na(filtered_df_sorted$gene_name))) {
    print("WARNING: Gene IDs missing - not provided by PRIDE")
    names(filtered_df_sorted)[names(filtered_df_sorted) == "gene_id"] <- "temp"
    names(filtered_df_sorted)[names(filtered_df_sorted) == "gene_name"] <- "gene_id"
    names(filtered_df_sorted)[names(filtered_df_sorted) == "temp"] <- "gene_name"
    
}

####################################
# save table
####################################

write.table(x = filtered_df_sorted, file = args[4], append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

