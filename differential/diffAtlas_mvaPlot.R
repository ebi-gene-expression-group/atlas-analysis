#!/usr/bin/env Rscript

# diffAtlas_mvaPlot()
#  - Create an "MvA" plot (average intensity against log2 fold change), and write it to a file.
# ARGUMENTS:
# 	- plotDataFile <- filename of differential expression stats to be used for plotting
# 	- contrastNAme <- text description of the contrast
# 	- plotFile <- filename to write the plot to
# 	- techType <- "microarray" or "rnaseq"

diffAtlas_mvaPlot <<- function(plotDataFile, contrastName, plotFile, techType) {

	# Read data
	madf <- read.delim(plotDataFile, stringsAsFactors=FALSE)

	# First filter out rows with any "NA"s in, can't use them.
	# complete.cases() function returns TRUE for rows where there is no missing
	# data (e.g. "NA"s).
	madf <- madf[complete.cases(madf),]

	requiredColumns <- c("avgExpr", "logFC", "adjPval")
	missingColumns <- setdiff(requiredColumns, colnames(madf))
	if(length(missingColumns) > 0) {
		stop(paste("Missing required MvA plot columns:", paste(missingColumns, collapse=", ")))
	}

	if(techType == "rnaseq") {
		# The RNA-seq plot uses a log10 x-axis, so zero-count rows cannot be plotted.
		madf <- madf[madf$avgExpr > 0,]
	}

	if(nrow(madf) == 0) {
		stop("No complete rows with plottable MvA data found.")
	}

	fdrCutoff <- 0.05
	foldChangeGuide <- 1

	madf$deCall <- "Non-DE"
	madf$deCall[madf$adjPval < fdrCutoff & madf$logFC > 0] <- "Up-regulated"
	madf$deCall[madf$adjPval < fdrCutoff & madf$logFC < 0] <- "Down-regulated"
	madf$deCall <- factor(madf$deCall, levels=c("Down-regulated", "Non-DE", "Up-regulated"))

	upCount <- sum(madf$deCall == "Up-regulated")
	downCount <- sum(madf$deCall == "Down-regulated")
	nonDeCount <- sum(madf$deCall == "Non-DE")
	summaryText <- paste0(
		"FDR < ", fdrCutoff, ": ",
		upCount, " up, ",
		downCount, " down, ",
		nonDeCount, " non-DE genes"
	)

	if(techType == "microarray") {

		# label for x-axis
		xAxisLabel = "average intensity"
	
	} else if(techType == "rnaseq") {
		
		# x-axis label different for RNA-seq data
		xAxisLabel = "average normalized count"
	}

	# load ggplot2
	library(ggplot2)

	nonDeData <- madf[madf$deCall == "Non-DE",]
	deData <- madf[madf$deCall != "Non-DE",]

	mvaPlot <- ggplot(madf, aes(x=avgExpr, y=logFC)) +
		geom_hline(yintercept=0, colour="grey35", linewidth=0.4) +
		geom_hline(yintercept=c(-foldChangeGuide, foldChangeGuide), colour="grey70", linetype="dotted", linewidth=0.35) +
		geom_point(data=nonDeData, aes(colour=deCall), alpha=0.28, size=0.55) +
		geom_point(data=deData, aes(colour=deCall), alpha=0.7, size=0.85) +
		scale_colour_manual(
			name=paste0("DE call (FDR < ", fdrCutoff, ")"),
			values=c("Down-regulated"="#2b6cb0", "Non-DE"="grey55", "Up-regulated"="#d73027"),
			drop=FALSE
		) +
		theme_minimal(base_size=12) +
		theme(
			panel.grid.minor=element_blank(),
			axis.line=element_line(colour="grey25", linewidth=0.35),
			legend.position="bottom",
			legend.direction="horizontal",
			legend.title=element_text(face="bold"),
			plot.title=element_text(face="bold", hjust=0.5, size=14),
			plot.subtitle=element_text(hjust=0.5, colour="grey35", margin=margin(b=10)),
			plot.caption=element_text(colour="grey40", size=9),
			plot.margin=margin(14, 18, 12, 14)
		)

	if(techType == "microarray") {
			
		mvaPlot <- mvaPlot +
			scale_y_continuous(breaks=pretty(madf$logFC, n=7)) +
			scale_x_continuous(breaks=pretty(madf$avgExpr, n=6))
	
	} else if(techType == "rnaseq") {
		
		# use log scale for DESeq results from RNA-seq data
		mvaPlot <- mvaPlot + scale_x_log10()
	} 
	
	# Label the axes. expression() function lets you make the "2" subscript (highly important! :)
	mvaPlot <- mvaPlot + xlab(xAxisLabel) + ylab(expression(log[2](fold~change))) +

		# Add the title (the description for this contrast). Use strwrap()
		# function to make it wrap after 70 characters.
		labs(
			title=paste(strwrap(contrastName, width=70), collapse="\n"),
			subtitle=summaryText,
			caption=paste0("Horizontal guide lines mark log2 fold-change 0 and +/-", foldChangeGuide, ".")
		)
	
	ggsave(filename=plotFile, plot=mvaPlot, width=7, height=6, dpi=300, bg="white")
}


# Run with arguments if there are any, otherwise don't do anything.

args <- commandArgs(TRUE)
if(length(args) > 0) {
	do.call(diffAtlas_mvaPlot, as.list(args))
}
