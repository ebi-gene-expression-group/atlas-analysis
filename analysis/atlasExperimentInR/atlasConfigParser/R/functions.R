# parseAtlasXML
# 	- Read Atlas XML config file and return a list of Analytics objects, as
# 	well as the experiment type.
# 	- The function returns a list with two elements: the list of Analytics
# 	objects, and the experiment type from the XML.
parseAtlasXML <- function( atlasXMLconfigFile ) {
	
	# Read the XML file.
	xmlTree <- xmlInternalTreeParse( atlasXMLconfigFile )

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

	# Since we can't return more than one thing, make a list to put the
	# analytics list in as well as the experiment type.
	parsedXMLlist <- list()
	
	parsedXMLlist$allAnalytics <- allAnalytics
	parsedXMLlist$experimentType <- atlasExperimentType

	return( parsedXMLlist )
}
