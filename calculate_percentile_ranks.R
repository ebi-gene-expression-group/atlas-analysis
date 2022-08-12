#!/usr/bin/env Rscript

# Load ExpressionAtlasInternal package.
suppressMessages( library( ExpressionAtlasInternal ) )

# Functions

# get_filename_start - Return the first part of the analytics file name --
# experiment accession and array design if appropriate.
get_filename_start <- function(analytics_filename) {
    analytics_basename <- basename(analytics_filename)
    filename_start <- sub("-analytics.tsv.unrounded", "", analytics_basename)
    return(filename_start)
}


# Script start

# Get commandline arguments.
args <- commandArgs(TRUE)

if (length(args) != 1) {
  stop("\nUsage:\n\tcalculate_percentile_ranks.R <analytics filename>\n\n")
}

# Get the analytics filename from the arguments.
analytics_filename <- args[1]
if (!grepl("-analytics.tsv.unrounded$", analytics_filename)) {
  stop(
    paste(analytics_filename,
    "does not look like a decorated and unrounded analytics filename. Please check."))
}

# Get the start of the filename (the experiment accession, and the array design
# accession if this is a microarray experiment).
filename_start <- get_filename_start(analytics_filename)

cat(paste("Reading analytics data from", analytics_filename, "...\n"))

# Read in analytics file.
analytics_data <- read.delim(analytics_filename,
                            header = TRUE, stringsAsFactors = FALSE)

cat("Successfully read analytics data.\n")

# Drop the Gene.Name column as we don't need it.
analytics_data$Gene.Name <- NULL
# Drop the Design.Element column if there is one as we don't need it.
if (any(grepl("Design.Element", colnames(analytics_data)))) {
  analytics_data$Design.Element <- NULL
}
# Also drop the p.value column as we don't need it either.
analytics_data <- analytics_data[,
                      which(!grepl("p.value", colnames(analytics_data)))]
# Also drop any t-statistic columns for the same reason.
if (any(grepl("t.statistic", colnames(analytics_data)))) {
  analytics_data <- analytics_data[,
                      which(!grepl("t.statistic", colnames(analytics_data)))]
}

# We need to write one percentile ranks file for each contrast.  Work out the
# contrast IDs.
stats_only_colnames <- colnames(analytics_data)[2:length(colnames(analytics_data))]
stats_only_colnames <- sub(".log2foldchange", "", stats_only_colnames)
contrast_ids <- unique(stats_only_colnames)

cat("Calculating percentile ranks for all contrasts...\n")

# Go through the contrast IDs and calculate the percentile ranks for the
# logFCs. This creates a list of data frames.
percentilesDataFrames <- lapply(contrast_ids, function(contrastID) {
  logFCcolname <- paste(contrastID, "log2foldchange", sep = ".")
  contrastDF <- data.frame(analytics_data$Gene.ID,
                           analytics_data[[logFCcolname]],
                           stringsAsFactors = FALSE)
  colnames(contrastDF) <- c("Gene.ID", logFCcolname)
  # Make sure we only take the best logFC for each gene, if there's more than
  # one entry per gene.  Sort the columns based on the gene ID and then
  # absolute logFC.
  contrastDF <- contrastDF[
                  order(contrastDF$Gene.ID, -abs(contrastDF[[logFCcolname]])), ]
  # Now just remove all the rows for gene IDs that come up as duplicates --
  # this leaves the first occurence of the gene ID (now always the largest
  # absolute logFC) and removes all the others.
  contrastDF <- contrastDF[which(!duplicated(contrastDF$Gene.ID)), ]
  # Now get the percentiles.
  percentiles <- quantile(abs(contrastDF[[logFCcolname]]),
                              probs = seq(0, 1, by = 0.01),
                              na.rm = TRUE)
  # Work out the indices of the percentiles vector that each logFC corresponds
  # to.
  percentileIndices <- findInterval(
                            abs(contrastDF[[logFCcolname]]),
                            sort(percentiles))

  # Get the numeric values from the names of the percentiles vector for each
  # logFC. These are the actual percentiles for each logFC.
  logFCpercentiles <- sub("%", "", names(percentiles)[percentileIndices])

  # Create a new data frame with the gene IDs and logFC percentiles.
  percentilesDF <- data.frame(Gene.ID = contrastDF$Gene.ID,
                              percentile = logFCpercentiles,
                              stringsAsFactors = FALSE)
  # Use the gene IDs as the row names.
  rownames(percentilesDF) <- percentilesDF$Gene.ID
  # Remove the gene IDs column.
  percentilesDF$Gene.ID <- NULL

  # Change the column name of the percentiles to include the contrast ID
  percentilesColname <- contrastID
  colnames(percentilesDF) <- percentilesColname

  # Sort the data frame rows based on the gene IDs.
  percentilesDF <- percentilesDF[order(rownames(percentilesDF)), , drop = FALSE]
  percentilesDF
})

cat("Successfully calculated percentile ranks.\n")

# Last step is to combine all the percentiles data frames into one, and write
# it out to a file.  First check that the gene IDs of the data frames are all
# the same.
geneIDref <- rownames(percentilesDataFrames[[1]])

invisible(sapply(percentilesDataFrames, function(percentilesDF) {
  if (all.equal(geneIDref, rownames(percentilesDF)) != TRUE) {
    stop("ERROR - Gene IDs of percentiles data frames are not identical, cannot join them together for writing.")
  }
}))

# If we're OK, join the data frames.
allPercentiles <- do.call("cbind", percentilesDataFrames)

# Add the gene IDs column before writing.
allPercentiles <- data.frame(Gene.ID = geneIDref, allPercentiles, stringsAsFactors = FALSE)

# Make the filename to write to.  FIXME: add the processing directory here.
percentilesFilename <- paste(filename_start, "-percentile-ranks.tsv", sep = "")

cat(paste("Writing percentile ranks to", percentilesFilename, "...\n"))

# Write the data frame.
write.table(allPercentiles, file = percentilesFilename, quote = FALSE, sep = "\t",
  row.names = FALSE)

cat("Successfully written percentile ranks.\n")