# getAtlasExperiment
#   - Download and return the SimpleList object representing a single
#   Expression Atlas experiment.
getAtlasExperiment <- function( experimentAccession ) {
    
    # Make sure we were given an experiment accession.
    if( missing( experimentAccession ) ) {

        stop( "Please provide experiment accession." )
    }

    # Make sure the experiment accession is in the correct format.
    if( !grepl( "^E-\\w{4}-\\d+$", experimentAccession ) ) {
        stop( paste( "\"", experimentAccession, "\" does not look like an ArrayExpress experiment accession. Cannot continue.", sep="" ) )
    }

    urlBase <- "http://www.ebi.ac.uk/gxa/experiments"

    atlasExperimentSummaryFile <- paste( experimentAccession, "-atlasExperimentSummary.Rdata", sep = "" )

    fullUrl <- paste( urlBase, experimentAccession, atlasExperimentSummaryFile, sep = "/" )
    
    cat( paste( "Loading Expression Atlas experiment summary from", fullUrl, "...\n" ) )

    result <- try( load( url( fullUrl ) ) )
    
    if( class( result ) == "try-error" ) {

        msg <- geterrmessage()

        stop( paste( "Error encountered while trying to download experiment summary:", msg ) )
    }

    # If we're still here, things must have worked ok.
    cat( paste( "Successfully downloaded experimentSummary object for", experimentAccession, "\n" ) )

    return( experimentSummary )
}
