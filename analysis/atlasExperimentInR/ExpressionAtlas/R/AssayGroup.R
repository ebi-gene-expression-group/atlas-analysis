# AssayGroup class
setClass( "AssayGroup", slots = c( assay_names = "vector", assay_group_id = "character", assay_group_label = "character" ) )


#####################
# AssayGroup Generics

# assay_group_id getter
setGeneric( "assay_group_id", function( object ) standardGeneric( "assay_group_id" ) )

# assay_names getter
setGeneric( "assay_names", function( object ) standardGeneric( "assay_names" ) )

# assay_group_label getter.
setGeneric( "assay_group_label", function( object ) standardGeneric( "assay_group_label" ) )

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
	
	# Get the assay_group element attributes.
	assayGroupAttrs <- xmlAttrs( assayGroupNode )
	
	# Get the assay group ID.
	assayGroupID <- assayGroupAttrs[[ "id" ]]

	# Get the assay group label
	assayGroupLabel <- tryCatch( 
		{
			assayGroupAttrs[[ "label" ]]
		},
		error = function( cond ) {
			message( paste( "Assay group", assayGroupID, "does not have a label" ) )
			return( NULL )
		}
	)

	# Add everything to the new object.
	.Object@assay_names <- assayGroupAssays
	.Object@assay_group_id <- assayGroupID

	if( !is.null( assayGroupLabel ) ) { 
		.Object@assay_group_label <- assayGroupLabel
	}

	return( .Object )
})

# assay_group_id getter
setMethod( "assay_group_id", "AssayGroup", function( object ) object@assay_group_id )

# Method for assay_names getter
setMethod( "assay_names", "AssayGroup", function( object ) object@assay_names )

# Method for assay_group_label getter
setMethod( "assay_group_label", "AssayGroup", function( object ) object@assay_group_label )
