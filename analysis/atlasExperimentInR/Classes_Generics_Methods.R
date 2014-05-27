#########
# Classes
#########

# AssayGroup
setClass( "AssayGroup", slots = c( assays = "vector", assay_group_id = "character" ) )








##########
# Generics
##########

# assay_group_id getter
setGeneric( "assay_group_id", function( object ) standardGeneric( "assay_group_id" ) )

# assay_group_id setter
setGeneric( "assay_group_id<-", function( object, value ) standardGeneric( "assay_group_id<-" ) )








#########
# Methods
#########

# assay_group_id getter
setMethod( "assay_group_id", "AssayGroup", function( object ) object@assay_group_id )

# assay_group_id setter
setReplaceMethod( "assay_group_id", 
	signature( 
		object="AssayGroup", value="character" 
	),
	function(object, value) {
		object@assay_group_id <- value
		return(object)
	}
)
	
