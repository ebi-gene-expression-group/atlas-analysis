# Contrast class
setClass( 
		 "Contrast", 
		 slots = c( 
				   contrast_id = "character", 
				   contrast_name = "character", 
				   reference_assay_group_id = "character", 
				   test_assay_group_id = "character",
				   batch_effects = "list"
) )

###################
# Contrast Generics

# contrast_id getter
setGeneric( "contrast_id", function( object ) standardGeneric( "contrast_id" ) )

# contrast_name getter
setGeneric( "contrast_name", function( object ) standardGeneric( "contrast_name" ) )

# reference_assay_group_id getter
setGeneric( "reference_assay_group_id", function( object ) standardGeneric( "reference_assay_group_id" ) )

# test_assay_group_id getter
setGeneric( "test_assay_group_id", function( object ) standardGeneric( "test_assay_group_id" ) )

# batch_effects getter
setGeneric( "batch_effects", function( object ) standardGeneric( "batch_effects" ) )

##################
# Contrast Methods

# constructor
setMethod( "initialize", "Contrast", function( .Object, contrastNode ) {

	# Get the contrast ID.
	contrastAttrs <- xmlAttrs( contrastNode )
	contrastID <- contrastAttrs[[ "id" ]]

	# Get the contrast name.
	contrastNameNode <- xmlElementsByTagName( contrastNode, "name" )$name
	contrastName <- xmlValue( contrastNameNode )

	# Get the reference and test assay group IDs.
	referenceAssayGroupIdNode <- xmlElementsByTagName( contrastNode, "reference_assay_group" )$reference_assay_group
	testAssayGroupIdNode <- xmlElementsByTagName( contrastNode, "test_assay_group" )$test_assay_group
	referenceAssayGroupId <- xmlValue( referenceAssayGroupIdNode )
	testAssayGroupId <- xmlValue( testAssayGroupIdNode )

	# Try for batch effects -- if there are any then create the object(s) and
	# add them as well.
	batchEffectObjects <- .create_batch_effect_objects( contrastNode )
	
	.Object@contrast_id <- contrastID
	.Object@contrast_name <- contrastName
	.Object@reference_assay_group_id <- referenceAssayGroupId
	.Object@test_assay_group_id <- testAssayGroupId	
	# Add batch effects if any.
	if( !is.null( batchEffectObjects ) ) {
		.Object@batch_effects <- batchEffectObjects
	}

	return( .Object )
})

# contrast_id getter
setMethod( "contrast_id", "Contrast", function( object ) object@contrast_id )

# contrast_name getter
setMethod( "contrast_name", "Contrast", function( object ) object@contrast_name )

# reference_assay_group_id getter
setMethod( "reference_assay_group_id", "Contrast", function( object ) object@reference_assay_group_id )

# test_assay_group_id getter
setMethod( "test_assay_group_id", "Contrast", function( object ) object@test_assay_group_id )

# batch_effects getter
setMethod( "batch_effects", "Contrast", function( object ) object@batch_effects )


# Look at the contrast node and create BatchEffect objects if any batch
# effects are found.
.create_batch_effect_objects <- function( contrastNode ) {
	
	# Get the batch_effects node, if it exists.
	batchEffectsNode <- xmlElementsByTagName( contrastNode, "batch_effects" )$batch_effects
	
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
