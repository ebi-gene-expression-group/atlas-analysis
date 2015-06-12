# Analytics class
setClass( 
		 "Analytics", 
		 slots = c( 
				   platform = "character", 
				   assay_groups = "list", 
				   atlas_contrasts = "list", 
				   batch_effects = "list" 
) )


####################
# Analytics Generics

# platform getter
setGeneric( "platform", function( object ) standardGeneric( "platform" ) )

# assay_groups getter
setGeneric( "assay_groups", function( object ) standardGeneric( "assay_groups" ) )

# atlas_contrasts getter
setGeneric( "atlas_contrasts", function( object ) standardGeneric( "atlas_contrasts" ) )

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
		platform <- gsub( "[ \t\n]", "", platform )
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

	# Get the atlas_contrasts node.
	contrastsList <- xmlElementsByTagName( analyticsNode, "contrasts" )
	contrastsNode <- contrastsList$contrasts

	# Get all the contrasts.
	allContrasts <- xmlElementsByTagName( contrastsNode, "contrast" )

	# Make a list containing Contrast objects.
	# Go through the contrast nodes.
	contrastObjects <- lapply( allContrasts, function( contrastNode ) {

		# Create a new contrast object.
		contrastObject <- new( "Contrast", contrastNode )
	})

	names( contrastObjects ) <- sapply( contrastObjects,
		function( contrast ) {
			contrast_id( contrast )
		}
	)

	# Try for batch effects -- if there are any then create the object(s) and
	# add them as well.
	batchEffectObjects <- .create_batch_effect_objects( analyticsNode )

	# Add everything we found to the object.
	.Object@platform <- platform
	.Object@assay_groups <- assayGroupObjects
	.Object@atlas_contrasts <- contrastObjects
	
	# Add batch effects if any.
	if( !is.null( batchEffectObjects ) ) {
		.Object@batch_effects <- batchEffectObjects
	}

	return( .Object )
})


# platform getter
setMethod( "platform", "Analytics", function( object ) object@platform )

# assay_groups getter
setMethod( "assay_groups", "Analytics", function( object ) object@assay_groups )

# atlas_contrasts getter
setMethod( "atlas_contrasts", "Analytics", function( object ) object@atlas_contrasts )

# getAllAssays
setMethod( "getAllAssays", "Analytics", function( object ) {
	
	assayGroups <- assay_groups( object )

	assay_names <- as.vector(
		unlist(
			sapply( assayGroups, 
				function( assayGroup ) {
					assay_names( assayGroup )
				}
			)
		)
	)

	names( assay_names ) <- NULL

	return( assay_names )
})


# Look at the analytics node and create BatchEffect objects if any batch
# effects are found.
.create_batch_effect_objects <- function( analyticsNode ) {
	
	# Get the batch_effects node, if it exists.
	batchEffectsNode <- xmlElementsByTagName( analyticsNode, "batch_effects" )$batch_effects
	
	# If a batch_effects node was found, make objects and return them.
	if( !is.null( batchEffectsNode ) ) {

		allBatchEffectNodes <- xmlElementsByTagName( batchEffectsNode, "batch_effect" )

		batchEffectObjects <- lapply( allBatchEffectNodes, 
			function( batchEffectNode ) {
				batchEffectObject <- new( "BatchEffect", batchEffectNode )
			}
		)

		return( batchEffectObjects )

	} else {
		
		# If no batch_effects object was found, just return null.
		return( NULL )
	}
}
