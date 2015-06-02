# BiologicalReplicate class.
setClass( "BiologicalReplicate", slots = c( assay_names = "vector", technical_replicate_id = "character" ) )

##############################
# BiologicalReplicate Generics

# assay_names getter.
setGeneric( "assay_names", function( object ) standardGeneric( "assay_names" ) )

# technical_replicate_id getter.
setGeneric( "technical_replicate_id", function( object ) standardGeneric( "technical_replicate_id" ) )

#############################
# BiologicalReplicate Methods

# Constructor.
setMethod( "initialize", "BiologicalReplicate", function( .Object, assayNames, techRepId = NULL ) {
	
	# Make sure that there's only one technical replicate ID.
	if( length( techRepId ) > 1 ) {
		cat( "ERROR - These do not look like technical replicate IDs:\n" )
		print( techRepId )
		cat( "ERROR - Please check technical replicate IDs in XML config.\n" )

		stop( "ERROR - Problem with technical replicate IDs. Cannot continue." )
	}

	.Object@assay_names <- assayNames
	.Object@technical_replicate_id <- techRepId

	return( .Object )
})


# Method for assay_names getter.
setMethod( "assay_names", "BiologicalReplicate", function( object ) object@assay_names )

# Method for technical_replicate_id getter.
setMethod( "technical_replicate_id", "BiologicalReplicate", function( object ) object@technical_replicate_id )
