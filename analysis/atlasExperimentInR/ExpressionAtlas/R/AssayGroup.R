# AssayGroup class
setClass( "AssayGroup", slots = c( biological_replicates = "vector", assay_group_id = "character", assay_group_label = "character" ) )


#####################
# AssayGroup Generics

# assay_group_id getter
setGeneric( "assay_group_id", function( object ) standardGeneric( "assay_group_id" ) )

# biological_replicates getter
setGeneric( "biological_replicates", function( object ) standardGeneric( "biological_replicates" ) )

# assay_names getter
setGeneric( "assay_names", function( object ) standardGeneric( "assay_names" ) )

# assay_group_label getter.
setGeneric( "assay_group_label", function( object ) standardGeneric( "assay_group_label" ) )

####################
# AssayGroup Methods

# constructor
setMethod( "initialize", "AssayGroup", function( .Object, assayGroupNode ) {

	# Get all the assay nodes.

	assayGroupAssays <- xmlElementsByTagName( assayGroupNode, "assay" )

	# Get a vector of biological replicates.
	biologicalReplicates <- .create_biological_replicates( assayGroupAssays )

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
			cat( paste( "Assay group", assayGroupID, "does not have a label\n" ) )
			return( NULL )
		}
	)

	# Add everything to the new object.
	.Object@biological_replicates <- biologicalReplicates
	.Object@assay_group_id <- assayGroupID

	if( !is.null( assayGroupLabel ) ) { 
		.Object@assay_group_label <- assayGroupLabel
	}

	return( .Object )
})

# Method for biological_replicates getter
setMethod( "biological_replicates", "AssayGroup", function( object ) object@biological_replicates )

# Method for assay_group_id getter
setMethod( "assay_group_id", "AssayGroup", function( object ) object@assay_group_id )

# Method for assay_names getter #FIXME
setMethod( "assay_names", "AssayGroup", function( object ) object@assay_names )

# Method for assay_group_label getter
setMethod( "assay_group_label", "AssayGroup", function( object ) object@assay_group_label )



######################
# Additional functions


.create_biological_replicates <- function( assayGroupAssays ) {
	
	# Initialise some variables needed later.
	bioRepsWithTechReps <- character()
	bioRepsWithoutTechReps <- character()

	# Get a vector of techincal replicate IDs, if there are any.
	technicalReplicateIds <- unique( unlist( sapply( assayGroupAssays, function( assayNode ) {
													assayNodeAttrs <- xmlAttrs( assayNode )
													assayNodeAttrs[[ "technical_replicate_id" ]]
												}
	) ) )
	
	# If we found some technical replicates, make BiologicalReplicates for
	# them.
	if( length( technicalReplicateIds ) > 0 ) {
		
		# First get the assay names for each technical replicate group.
		# Use the "simplify = FALSE" setting in sapply() so that the named list
		# structure is not lost.
		techRepAssays <- sapply( technicalReplicateIds, simplify = FALSE, function( techRepId ) {
								techRepAssays <- .find_assays_for_tech_rep_id( allAssays, techRepId )
		})
		
		# Now make a vector of BiologicalReplicate objects.
		bioRepsWithTechReps <- unlist( sapply( names( techRepAssays ), function( techRepId ) {
			
			assayNames <- techRepAssays[[ techRepId ]]

			biologicalReplicate <- new( "BiologicalReplicate", assayNames, techRepId )
	
		} ) )
	}

	# Get the names of the assays that are not technical replicates.
	noTechRepAssays <- unlist( sapply( assayGroupAssays, function( assayNode ) {
									  assayNodeAttrs <- xmlAttrs( assayNode )
									  if( length( assayNodeAttrs ) == 0 ) { xmlValue( assayNode ) }
									}
	) )

	# Now create the BiologicalReplicate objects for the assays that are not
	# technical replicates.
	bioRepsWithoutTechReps <- sapply( noTechRepAssays, function( assayName ) {
								biologicalReplicate <- new( "BiologicalReplicate", assayName )
	})
	
	# Return the BiologicalReplicate objects that were created.
	return( c( bioRepsWithTechReps, bioRepsWithoutTechReps ) )
}


.find_assays_for_tech_rep_id <- function( assayGroupAssays, tech_rep_id ) {

	techRepAssays <- sapply( assayGroupAssays, function( assayNode ) {
		assayNodeAttrs <- xmlAttrs( assayNode )
		if( length( assayNodeAttrs ) > 0 ) {
			if( assayNodeAttrs[[ "technical_replicate_id" ]] == tech_rep_id ) { xmlValue( assayNode ) }
		}
	})

	return( unlist( techRepAssays ) )
}
