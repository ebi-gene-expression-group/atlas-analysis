# AssayGroup class
setClass( "AssayGroup", slots = c( assay_names = "vector", assay_group_id = "character" ) )


#####################
# AssayGroup Generics

# assay_group_id getter
setGeneric( "assay_group_id", function( object ) standardGeneric( "assay_group_id" ) )

# assay_names getter
setGeneric( "assay_names", function( object ) standardGeneric( "assay_names" ) )

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
	
	# Make assay names "R-safe" as they have to match ones in column headings
	# later on.
	assayGroupAssays <- make.names( assayGroupAssays )

	# Remove the names from the vector of assay names.
	names(assayGroupAssays) <- NULL

	# Get the assay group ID.
	assayGroupAttrs <- xmlAttrs( assayGroupNode )
	assayGroupID <- assayGroupAttrs[[ "id" ]]

	.Object@assay_names <- assayGroupAssays
	.Object@assay_group_id <- assayGroupID
	return( .Object )
})

# assay_group_id getter
setMethod( "assay_group_id", "AssayGroup", function( object ) object@assay_group_id )

# assay_names
setMethod( "assay_names", "AssayGroup", function( object ) object@assay_names )
