\name{summarizeAtlasExperiment}
\alias{summarizeAtlasExperiment}
\title{Summarize an Expression Atlas experiment with Bioconductor objects}
\description{
	Create a \code{SimpleList} object containing \code{ExpressionSet},
	\code{SummarizedExperiment}, and/or \code{MAList} objects storing data from
	an Expression Atlas experiment.
}
\usage{
	summarizeAtlasExperiment( experimentAccession, atlasExperimentDirectory )
}
\arguments{
	\item{experimentAccession}{
	A character string containing the ArrayExpress accession of the experiment
	}
	\item{atlasExperimentDirectory}{
	Path to the directory containing the Atlas XML config file and expressions
	matrices
	}
}
\value{
	A SimpleList object representing one Expression Atlas experiment. Each
	element is either an ExpressionSet, SummarizedExperiment, or MAList object.

	For microarray experiments, the elements are named according to their
	ArrayExpress array design accession. For RNA-seq experiments, there will
	only be one element, named "rnaseq".

	SummarizedExperiment objects are used for RNA-seq data. They contain
	un-normalized, raw counts created with HTSeq. ExpressionSet objects are
	used for 1-colour microarray data. MAList objects are used for 2-colour
	microarray data. Microarray objects contain normalized intensity values.
	Affymetrix intensities are normalized using RMA, from the oligo package.
	Agilent 1-colour array intensities are normalized using between-array
	quantile normalization in the limma package. Agilent 2-colour array
	intensities are normalized using within-array loess normalization in limma.
}
\examples{

# Summarize an Atlas microarray experiment.
arrayExperimentSummaryList <- summarizeAtlasExperiment( 
	"E-MEXP-1025", 
	file.path( "extdata", "E-MEXP-1025" ) 
)

# Get some ExpressionSet objects.
expressionSet_affy18 <- arrayExperimentSummaryList[[ "A-AFFY-18" ]]
expressionSet_affy19 <- arrayExperimentSummaryList[[ "A-AFFY-19" ]]

# See the normalization method used.
preproc( experimentData( expressionSet_affy18 ) )

# Summarize an Atlas RNA-seq experiment.
rnaseqExperimentSummaryList <- summarizeAtlasExperiment(
	"E-GEOD-38400",
	file.path( "extdata", "E-GEOD-38400" )
)

# Get the SummarizedExperiment object.
summarizedExperiment <- rnaseqExperimentSummaryList$rnaseq

# Get the analysis methods.
analysisMethods <- exptData( summarizedExperiment )

# See what methods are documented.
names( analysisMethods )

# See the filtering steps used.
analysisMethods$filtering
}
