#!/usr/bin/env Rscript
# 
# Calculate Pearson's product-moment correlation coefficient or Spearman Rank
# correlation using two sets of RNA-seq read counts (or *PKMs) produced by
# different quantification methods.  Produce a scatterplot and display
# correlation coefficient on it.
#
# Provide two tab-delimited text files with genes in rows and sequencing runs
# (or averaged counts/FPKMs for a group of replicates) in columns.
# Provide gene IDs in a file and/or column headings to subset from the above
# files. To use all genes and/or all runs, use "all" instead (see below
# examples).
# 
# For each tab-delimited file, counts(/*PKMs) from all columns specified are
# summed for each gene. The resulting two vectors of counts are correlated and
# plotted.
#
# Examples:
# 1. All genes and all columns:	
# 	Rscript correlate.R all <counts file 1> all <counts file 2> all <plot filename> pearson[|spearman]
#
# 2. Only some genes but all columns:
# 	Rscript correlate.R <gene IDs file> <counts file 1> all <counts file 2> all <plot filename> pearson[|spearman] 
#
# 3. All genes but only some columns:
# 	(a) sequencing runs from the same experiment.
# 		Rscript correlate.R all <counts file 1> SRR001,SRR002 <counts file 2> SRR001,SRR002 <plot filename> pearson[|spearman]
# 	(b) averaged biological replicates from different experiments (e.g. heart
# 	was in column "g8" in one experiment and "g10" in the other)
# 		Rscript correlate.R all <counts file 1> g8 <counts file 2> g10 <plot filename> pearson[|spearman]


# correlate
# ARGUMENTS:
# 	- genes : either "all" or name of a file containing a list of gene IDs (one per line) to use in correlation.
# 	- countsFile1 : name of first file containing data to use in correlation.
# 	- cols1 : comma-separated list of column headings to select from file 1.
# 	- countsFile2 : name of second file containing data to use in correlation.
# 	- cols2 : comma-separated list of column headings to select from file 2.
# 	- outFile : name of PDF file to write plot of correlation to.
# 	- corType : either "pearson" for Pearson's product-moment correlation, or "spearman" for Spearman Rank correlation.
correlate <<- function(genes, countsFile1, cols1, countsFile2, cols2, outFile, corType, makePlot = "plot") {
	
    if( makePlot != "plot" && makePlot != "noplot" ) {
        stop( "makePlot value must be either plot or noplot" )
    }

	# Get labels for plot axes from filenames.
	counts1Lab <- basename(countsFile1)
	counts2Lab <- basename(countsFile2)

	# Read in counts.
	counts1 <- read.delim(countsFile1, header=TRUE, stringsAsFactors=FALSE)
	rownames(counts1) <- counts1[,1]
	counts1[,1] <- NULL
	
	counts2 <- read.delim(countsFile2, header=TRUE, stringsAsFactors=FALSE)
	rownames(counts2) <- counts2[,1]
	counts2[,1] <- NULL

	# Subset selected columns from the first file if need be.
	if(cols1 != "all") {
		
		cols1 <- unlist(strsplit(cols1, ","))

		counts1 <- counts1[,cols1, drop=FALSE]
	} else {
		message(paste("Using all columns from", countsFile1, "\n"))
	}
	
	# Subset selected columns from the second file if need be.
	if(cols2 != "all") {
		
		cols2 <- unlist(strsplit(cols2, ","))

		counts2 <- counts2[,cols2, drop=FALSE]
	} else {
		message(paste("Using all columns from", countsFile2, "\n"))
	}

	# Read genes from a file unless "all" is passed.
	if(genes != "all") {
		
		message(paste("Using genes from", genes, "\n"))
		genes <- scan(file=genes, what="character")
		
		counts1 <- counts1[genes, , drop=FALSE]

		# Check for genes that weren't found in the first counts matrix.
		if(length(grep("NA", rownames(counts1))) > 0 ) {
			# Get IDs of genes that weren't found.
			unwantedGenes1 <- genes[grep("NA", rownames(counts1))]
			
			# Log that they weren't found and will be removed.
			message(paste("The following genes were not found in", countsFile1, "\n"))
			message(paste(unwantedGenes1, collapse="\n"))
			message("\nThey will be removed.\n\n")
			
			# Remove them from counts1.
			counts1 <- counts1[-grep("NA", rownames(counts1)),,drop=FALSE]

			# Also remove them from counts2 if necessary.
			if(length(which(rownames(counts2) %in% unwantedGenes1)) > 0) {
				counts2 <- counts2[-which(rownames(counts2) %in% unwantedGenes1),,drop=FALSE]
			}
		}
		
		counts2 <- counts2[genes, , drop=FALSE]
		
		# Check for genes that weren't found in the second counts matrix.
		if(length(grep("NA", rownames(counts2))) > 0 ) {
			# Get IDs of genes that weren't found.
			unwantedGenes2 <- genes[grep("NA", rownames(counts2))]
			
			# Log that they weren't found and will be removed.
			message(paste("The following genes were not found in", countsFile2, "\n"))
			message(paste(unwantedGenes2, collapse="\n"))
			message("\nThey will be removed.\n\n")
			
			# Remove them from counts2.
			counts2 <- counts2[-grep("NA", rownames(counts2)),,drop=FALSE]

			# Also remove them from counts1 if necessary.
			if(length(which(rownames(counts1) %in% unwantedGenes2)) > 0) {
				counts1 <- counts1[-which(rownames(counts1) %in% unwantedGenes2),,drop=FALSE]
			}
		}
			
	} else {
		message("\nUsing all genes\n")
	}

	# Sort rows by gene IDs.
	counts1 <- counts1[order(rownames(counts1)),,drop=FALSE]
	counts2 <- counts2[order(rownames(counts2)),,drop=FALSE]

	# Make sure that the row names in each set of counts are identical. It
	# doesn't make sense to try and correlate counts from one gene from one
	# method with counts from a totally different gene from the other method.
	stopifnot(all(rownames(counts1) == rownames(counts2)))
	
	# Add counts for each gene and put the two vectors into a data frame.
	counts1Sum <- apply(counts1, 1, function(x) { sum(x) })
	counts2Sum <- apply(counts2, 1, function(x) { sum(x) })
	sumDF <- data.frame(counts1Sum=counts1Sum, counts2Sum=counts2Sum)

	# Correlate and plot.
	cor.scatterplot(sumDF, "counts1Sum", "counts2Sum", counts1Lab, counts2Lab, outFile, corType, makePlot )
}


# Correlation plot function
# ARGUMENTS:
# 	- df : data frame containing sum of counts from the two matrices supplied.
# 	- colA : name of first column in above data frame.
# 	- colB : name of second column in above data frame.
# 	- counts1Lab : name of first data file.
# 	- counts2Lab : name of second data file.
# 	- outFile : name of PDF file to write to.
# 	- corType : "pearson" or "spearman".
cor.scatterplot <- function(df, colA, colB, counts1Lab, counts2Lab, outFile, corType, makePlot) {
  	
	# Check corType makes sense
	if(corType != "pearson" && corType != "spearman") {
		stop(paste("Unknown correlation type:", corType))
	}

	# add 0.1 to all counts (why??).
	del <- 0.1

	message(paste("Running correlation with method =", corType, "\n"))

	# use cor() to get Pearson's or Spearman's correlation coefficient, and round it to 2 decimal places.
    x<-round(cor(df[,colA], df[,colB], method=corType), 2)

    if( makePlot == "noplot" ) {
        cat( paste( x, "\n" ) )
    }
    else {
	
        # Start PDF device to write plot to.
        pdf(file=outFile)

        # Plot!
        # Name for plot depends on correlation type:
        if(corType == "pearson") {
            plotNameStart <- "Pearson's coefficient = "
        } else {
            plotNameStart <- "Spearman's coefficient = "
        }
        plot(df[,colA]+del, df[,colB]+del, log="xy", type="p", pch=".", xlab=counts1Lab, ylab=counts2Lab, main=paste(plotNameStart, x, sep=""))
        # add line y=x
        abline(0, 1, col="red")
        # close PDF.
        dev.off()
    }
}



# Call correlate() function if we've got arguments.
args <- commandArgs(TRUE)
if(length(args) > 0) {
	do.call(correlate, as.list(args))
}
