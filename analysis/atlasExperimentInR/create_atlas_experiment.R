
source("/ebi/microarray/home/mkeays/Atlas/git/atlasprod/analysis/atlasExperimentInR/AssayGroup.R")
source("/ebi/microarray/home/mkeays/Atlas/git/atlasprod/analysis/atlasExperimentInR/Analytics.R")


parseAtlasXML <- function(filename) {
	
	# Load XML package.
	library( XML )

	# Read the XML file.
	xmlTree <- xmlInternalTreeParse( filename )

	# Get the configuration node -- the root of the tree.
	configNode <- xmlRoot( xmlTree )

	# Get the config node attributes.
	configNodeAttrs <- xmlAttrs( configNode )

	# Get the Atlas experiment type. We'll use this later to decide whether to
	# create an ExpressionSet/MAList or SummarizedExperiment.
	atlasExperimentType <- configNodeAttrs[[ "experimentType" ]]
	
	# Go through the analytics node(s).
	# Get them from the configuration node.
	allAnalyticsNodes <- xmlElementsByTagName( configNode, "analytics" )

	# Create a list of Analytics objects.
	allAnalytics <- lapply( allAnalyticsNodes, function( analyticsNode ) {
		analyticsObject <- new( "Analytics", atlasExperimentType, analyticsNode )
	})
	
	names( allAnalytics ) <- sapply( allAnalytics,
		function( analytics ) {
			platform( analytics )
		}
	)

	return( allAnalytics )
}


