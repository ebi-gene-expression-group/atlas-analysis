# BatchEffect class
setClass( "BatchEffect", slots = c( batches = "list", effect_name = "character" ) )

######################
# BatchEffect Generics

# batches getter
setGeneric( "batches", function( object ) standardGeneric( "batches" ) )

# effect_name getter
setGeneric( "effect_name", function( object ) standardGeneric( "effect_name" ) )

#####################
# BatchEffect Methods

# constructor
setMethod( "initialize", "BatchEffect", function( .Object, batchEffectNode ) {
	
	# Get the batch name.
	batchName <- xmlAttrs( batchEffectNode )[[ "name" ]]

	# Get the batch nodes.
	batchNodes <- xmlElementsByTagName( batchEffectNode, "batch" )
	
	# Go through the batches and make a list indexing the assay names by their
	# batch value. The list should end up like e.g.
	#
	# allBatches$male = [ assay1, assay2, assay3, ... ]
	# allBatches$female = [ assay8, assay9, assay10, ... ]
	allBatches <- sapply( batchNodes,
		function( batchNode ) {
			
			# Get the value for this batch, e.g. male, female, ...
			batchValue <- xmlAttrs( batchNode )[[ "value" ]]
			
			# Get all the assay nodes.
			allAssayNodes <- xmlElementsByTagName( batchNode, "assay" )
			
			# Get a vector of the assay names.
			assayNames <- sapply( allAssayNodes, 
				function( assayNode ) {
					make.names( xmlValue( assayNode ) )
				}
			)
			
			# Remove the names from the assay names vector as they're useless.
			names( assayNames ) <- NULL
			
			# Create a list for the batch.
			batchList <- list()
			
			# Add the assay names vector under this batch's value.
			batchList[[ batchValue ]] <- assayNames
			
			# Return the batch list.
			batchList
		}
	)

	# All the names in the finished list have "batch." added at the beginning
	# so remove this.
	names( allBatches ) <- sub( "^batch.", "", names( allBatches ) )
	
	# Add the batch name and batches list to the object and return it.
	.Object@effect_name <- batchName
	.Object@batches <- allBatches

	return( .Object )
})

# batches getter method
setMethod( "batches", "BatchEffect", function( object ) object@batches )

# effect_name getter method
setMethod( "effect_name", "BatchEffect", function( object ) object@effect_name )
