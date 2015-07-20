# Contrast class
setClass( "Contrast", slots = c( contrast_id = "character", contrast_name = "character", reference_assay_group_id = "character", test_assay_group_id = "character" ) )

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
	
	.Object@contrast_id <- contrastID
	.Object@contrast_name <- contrastName
	.Object@reference_assay_group_id <- referenceAssayGroupId
	.Object@test_assay_group_id <- testAssayGroupId

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

