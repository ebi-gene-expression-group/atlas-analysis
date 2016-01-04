\name{parseAtlasConfig}
\alias{parseAtlasConfig}
\title{Read XML config from an Expression Atlas experiment}
\description{
	Read an XML config from an Expression Atlas experiment and create a list
	containing information on the analytics and assay group elements detailed
	therein.
}
\usage{
	parseAtlasConfig( atlasXMLconfigFile )
}
\arguments{
	\item{ atlasXMLconfigFile }{
	Path to the Expression Atlas XML config file.
	}
}
\value{
	A list with two elements: the experiment type, and a list of all Analytics
	elements found in the XML config file.
}
\examples{
\dontrun{
	# Path to atlas XML config file.
	atlasXMLconfigFile <- system.file( "extdata", "E-MTAB-1066-configuration.xml", package = "ExpressionAtlasInternal" )

	# Parse the config file into a list
	parsedXMLlist <- parseAtlasConfig( atlasXMLconfigFile )

	# Get the experiment type
	atlasExperimentType <- parsedXMLlist$experimentType

	# Get the Analytics objects (which contain AssayGroups)
	allAnalytics <- parsedXMLlist$allAnalytics
}
}
