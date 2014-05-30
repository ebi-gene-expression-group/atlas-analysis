source("/ebi/microarray/home/mkeays/Atlas/git/atlasprod/analysis/atlasExperimentInR/AssayGroup.R")

# Analytics class
setClass( "Analytics", slots = c( platform = "character", assay_groups = "list" ) )


####################
# Analytics Generics

# platform getter
setGeneric( "platform", function( object ) standardGeneric( "platform" ) )

# assay_groups getter
setGeneric( "assay_groups", function( object ) standardGeneric( "assay_groups" ) )

# getAllAssays
setGeneric( "getAllAssays", function( object ) standardGeneric( "getAllAssays" ) )

###################
# Analytics Methods

# constructor
setMethod( "initialize", "Analytics", function( .Object, atlasExperimentType, analyticsNode ) {

	# If this is a microarray experiment, get the array design.
	if( grepl( "array", atlasExperimentType ) ) {

		# Get a list of array design nodes.
		arrayDesignNodeList <- xmlElementsByTagName( analyticsNode, "array_design" )

		# Get the array design node.
		arrayDesignNode <- arrayDesignNodeList$array_design

		# Get the array design accession.
		platform <- xmlValue( arrayDesignNode )

		# Sometimes the array design has extra tabs and newlines that we don't
		# want, so remove them.
		platform <- gsub( "[\t\n]", "", platform )
	}
	else {
		platform <- "rnaseq"
	}

	# Now get the assay_groups node.
	assayGroupsList <- xmlElementsByTagName( analyticsNode, "assay_groups" )
	assayGroupsNode <- assayGroupsList$assay_groups

	# Get all the assay groups
	allAssayGroups <- xmlElementsByTagName( assayGroupsNode, "assay_group" )

	# Make a list containing AssayGroup objects.
	# Go through the assay group nodes.
	assayGroupObjects <- lapply( allAssayGroups, function( assayGroupNode ) {

		# Create a new AssayGroup object.
		assayGroupObject <- new( "AssayGroup", assayGroupNode )
	})

	names( assayGroupObjects ) <- sapply( assayGroupObjects,
		function( assayGroup ) {
			assay_group_id( assayGroup )
		}
	)

	.Object@platform <- platform
	.Object@assay_groups <- assayGroupObjects
	return( .Object )
})


# platform getter
setMethod( "platform", "Analytics", function( object ) object@platform )

# assay_groups getter
setMethod( "assay_groups", "Analytics", function( object ) object@assay_groups )

# getAllAssays
setMethod( "getAllAssays", "Analytics", function( object ) {
	
	assayGroups <- assay_groups( object )

	assays <- as.vector(
		unlist(
			sapply( assayGroups, 
				function( assayGroup ) {
					assays( assayGroup )
				}
			)
		)
	)

	names( assays ) <- NULL

	return( assays )
})

