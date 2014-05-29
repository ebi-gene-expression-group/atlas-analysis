# AssayGroup class
setClass( "AssayGroup", slots = c( assays = "vector", assay_group_id = "character" ) )


#####################
# AssayGroup Generics

# assay_group_id getter
setGeneric( "assay_group_id", function( object ) standardGeneric( "assay_group_id" ) )

# assays getter
setGeneric( "assays", function( object ) standardGeneric( "assays" ) )

####################
# AssayGroup Methods

# constructor
setMethod( "initialize", "AssayGroup", function( .Object, assayGroupNode ) {

	# Get all the assay nodes.

	allAssays <- xmlElementsByTagName( assayGroupNode, "assay" )

	# Get a vector of the assay names.
	assayGroupAssays <- sapply( allAssays, function( assayNode ) {
		xmlValue( assayNode )
	})

	# Remove the names from the vector of assay names.
	names(assayGroupAssays) <- NULL

	# Get the assay group ID.
	assayGroupAttrs <- xmlAttrs( assayGroupNode )
	assayGroupID <- assayGroupAttrs[[ "id" ]]

	.Object@assays <- assayGroupAssays
	.Object@assay_group_id <- assayGroupID
	return( .Object )
})

# assay_group_id getter
setMethod( "assay_group_id", "AssayGroup", function( object ) object@assay_group_id )

# assays
setMethod( "assays", "AssayGroup", function( object ) object@assays )
