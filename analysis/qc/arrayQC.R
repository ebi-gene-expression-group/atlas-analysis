#!/ebi/microarray/home/biocep/local/lib64/R/bin/Rscript
#
# Run quality assessment on microarray data prior to loading to Atlas, using
# the arrayQualityMetrics Bioconductor package. This script requires a
# tab-delimited annotation file mapping factor values to raw data filenames,
# created by arrayQC.pl.

# atlasArrayQC
# 	- Run arrayQualityMetrics for an experiment-array design.
# Arguments:
# 	- annotationFile : filename of annotation file containing mapping between
# 	factor values and raw data files.
# 	- exptType : affy, agil1 or agil2
# 	- exptAcc : experiment accession
# 	- arrayDesign : array design accession
# 	- outDir : directory to write quality report into
# 	- miRBaseFile : path to file with miRBase probe name mappings (0 for a
# 	non-miRNA experiment)
atlasArrayQC <- function(annotationFile, exptType, exptAcc, arrayDesign, outDir, miRBaseFile) {
	# Make a title for the report
	reportTitle <- paste("Microarray quality metrics report for", exptAcc, "on array design", arrayDesign)
	
	# The name of the column with factor values in the annotations file. This
	# will be "FactorValue" for one-colour experiments but will change to "Cy5"
	# for two-colour ones.  The values in this column are used to colour points
	# in the plots in the report.
	factorColName <- "FactorValue"

	# Read in raw data and add factor values to object created.
	# Affymetrix data
	if(exptType == "affy") {
		# Read in the annotation file with factor values and corresponding raw data
		# files.
		annotations <- read.delim(annotationFile, header = TRUE, stringsAsFactors = FALSE)
		rownames(annotations) <- annotations$AssayName

		# Create Biobase AnnotatedDataFrame object with annotations.
		library(Biobase)
		annotDF <- new("AnnotatedDataFrame", annotations)
	
		# Load oligo package.
		library(oligo)
		
		# Read files into ExpressionFeatureSet object.
		dataSet <- try({read.celfiles(annotations$FileName)})
		if(class(dataSet) == "try-error") {
			return(dataSet)
		}
		
		# Add annotations to the object.
		phenoData(dataSet) <- annotDF
	}
	# Agilent data.
	# ***For miRNA do we need to subset probes prior to QM calc??
	# One-colour Agilent data
	else if(exptType == "agil1") {
		# Read in the annotation file with factor values and corresponding raw data
		# files.
		annotations <- read.delim(annotationFile, header = TRUE, stringsAsFactors = FALSE)
		rownames(annotations) <- annotations$AssayName

		# Create Biobase AnnotatedDataFrame object with annotations.
		library(Biobase)
		annotDF <- new("AnnotatedDataFrame", annotations)
	
		# Load limma package
		library(limma)
		
		# Read in the data to an "EListRaw" object.
		dataSet <- try({read.maimages(annotations$FileName, source="agilent", green.only=TRUE)})
		
		# If we have a miRBase mapping file, subset the data to only include
		# probes that are there.
		if(miRBaseFile != 0) {
			dataSet <- subsetProbes(dataSet, miRBaseFile)
		}
	
		# Make it into an "NChannelSet" object for arrayQualityMetrics to read.
		# For this need to put the expressions into an AssayData object.
		aData <- assayDataNew(storage.mode = "lockedEnvironment", exprs = dataSet$E)
		# Add the rownames from the annotations data frame (filenames) as the sample names.
		sampleNames(aData) <- rownames(annotations)
		# Create the NChannelSet object with the expressions and annotations.
		dataSet <- new("NChannelSet", assayData = aData, phenoData = annotDF)
		
		# Check things worked and return the error if not.
		if(class(dataSet) == "try-error") {
			return(dataSet)
		}
	}
	# Two-colour Agilent data
	else if(exptType == "agil2") {
		# Use the "Cy5" column from the annotations to colour points in the
		# report.
		factorColName <- "Cy5"

		# Load limma
		library(limma)

		# Create a "targets" data frame to show which samples are labelled with
		# which dye in the report.
		targets <- readTargets(annotationFile)
		
		# Read in the data to an "RGList" object.
		dataSet <- try({read.maimages(targets, source="agilent")})
		
		# Check things worked and return the error if not.
		if(class(dataSet) == "try-error") {
			return(dataSet)
		}
		
		# If we have a miRBase mapping file, subset the data to only include
		# probes that are there.
		if(miRBaseFile != 0) {
			dataSet <- subsetProbes(dataSet, miRBaseFile)
		}
	}
	# Something we don't recognise.
	else { stop(paste("Can't handle", exptType, "experiments yet", sep=" ")) }

	# Load arrayQualityMetrics package
	library(arrayQualityMetrics)
	
	# Run quality metrics calculation
	arrayQualityMetrics(expressionset = dataSet, 
						outdir = outDir,
						force = TRUE,
						do.logtransform = TRUE,
						intgroup = factorColName,
						spatial = FALSE,
						reporttitle = reportTitle)
	
}


# Copied this from normalizeOneExperiment.R but should prob take it out to a
# separate script to source.
# subsetProbes
#  - For microRNA, we only want to normalize using probesets that are mapped to
#  the latest release of miRBase (www.mirbase.org). Subset the data here.
# ARGUMENTS:
# 	- dataSet <- for now either EListRaw or RGList object.
# 	- miRBaseFile <- filename for miRBase mappings.
subsetProbes <<- function(dataSet, miRBaseFile) {

	print(paste("There are", nrow(dataSet), "rows of data"))
	
	print(paste("Subsetting data for probes found in", miRBaseFile))

	# Read file with miRBase mappings into a data frame.
	miRBaseProbeMapping <- read.delim(miRBaseFile, stringsAsFactors = FALSE)
	
	# Probesets can be repeated; take unique set from design_element column in
	# miRBase mapping data frame.
	design_element <- unique(miRBaseProbeMapping$design_element)

	# Subset the data, taking only probesets from miRBase mapping file.
	# Agilent data is read in to one of two classes:
	# 	- "EListRaw" : 1-colour
	# 	- "RGList" : 2-colour
	if(class(dataSet) %in% c("EListRaw", "RGList")) { 
		
		# Subset dataSet using design_element column from miRBase mapping file.
		dataSet <- dataSet[which(dataSet$genes$ProbeName %in% design_element), ]
		
		# Check if we have any rows left after subsetting.
		if(nrow(dataSet) == 0) {
			# If not, that means no probe names matched any of the values in
			# the design_element column, so die.
			stop("No probes in raw data match values in design_element column of miRBase mapping file. Cannot proceed.")
		}
	}
	# If it's not an "ElistRaw" or "RGList" object we can't handle it (yet).
	else { print("Don't know how to subset for this type of array") }

	print(paste("After subsetting there are", nrow(dataSet), "rows remaining"))

	# TODO: add something here to subset for Affymetrix arrays.
	# Need to think about how to map probes to probesets because before
	# normalization probeset info is not available?
	
	return(dataSet)
}


# Run with arguments if there are any, otherwise don't do anything.
args <- commandArgs(TRUE)
if(length(args) > 0) {
	do.call(atlasArrayQC, as.list(args))
}
