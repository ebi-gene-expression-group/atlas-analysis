# getAtlasExperiment
#   - Download and return the SimpleList object representing a single
#   Expression Atlas experiment.
getAtlasExperiment <- function( experimentAccession ) {
    
    # Make sure the experiment accession is in the correct format.
    if( ! .isValidExperimentAccession( experimentAccession ) ) {

        stop( "Experiment accession not valid. Cannot continue." )
    }

    urlBase <- "http://www.ebi.ac.uk/gxa/experiments"

    atlasExperimentSummaryFile <- paste( 
        experimentAccession,
        "-atlasExperimentSummary.Rdata", 
        sep = "" 
    )

    fullUrl <- paste( 
        urlBase, 
        experimentAccession, 
        atlasExperimentSummaryFile, 
        sep = "/" 
    )
    
    cat( 
        paste( 
            "Downloading Expression Atlas experiment summary from:\n", 
            fullUrl, 
            "...\n" 
        ) 
    )

    loadResult <- try( load( url( fullUrl ) ) )
    
    if( class( loadResult ) == "try-error" ) {

        msg <- geterrmessage()

        stop( 
            paste( 
                "Error encountered while trying to download experiment summary:", 
                msg 
            ) 
        )
    }
    
    # Make sure experiment summary object exists before trying to return it.
    getResult <- try( get( "experimentSummary" ) )

    if( class( getResult ) == "try-error" ) {
        
        stop( 
            "ERROR - Download appeared successful but no experiment summary object was found." 
        )
    }

    # If we're still here, things must have worked ok.
    cat( 
        paste( 
            "Successfully downloaded experiment summary object for", 
            experimentAccession, 
            "\n" 
        ) 
    )

    return( get( "experimentSummary" ) )
}


# getAtlasData
#   - Download SimpleList objects for one or more Expression Atlas experiments
#   and return then in a list.
getAtlasData <- function( experimentAccessions ) {

    if( missing( experimentAccessions ) ) {

        stop( "Please provide a vector of experiment accessions to download." )
    }

    # Make sure experimentAccessions is a vector.
    if( ! is.vector( experimentAccessions ) ) {
        
        stop( "Please provide experiment accessions as a vector." )
    }

    # Go through each one and download it, creating a list.
    # So that the list has the experiment accessions as the names, use them as
    # names for the vector before starting to create the list.
    names( experimentAccessions ) <- experimentAccessions

    experimentSummaryList <- lapply( experimentAccessions, function( experimentAccession ) {

        if( .isValidExperimentAccession( experimentAccession ) ) {
            
            getAtlasExperiment( experimentAccession )

        } else {

            warning( 
                paste( 
                    "Not attempting to download data for invalid accession \"", 
                    experimentAccession, 
                    "\"", 
                    sep="" 
                ) 
            )
        }
    } )

    # Make sure we got something back -- if all the accessions were invalid
    # then experimentSummaryList may be empty.

    return( experimentSummaryList )
}


# .isValidExperimentAccession
#   - Return TRUE if experiment accession matches expected ArrayExpress
#   experiment accession pattern. Return FALSE otherwise.
.isValidExperimentAccession <- function( experimentAccession ) {

    if( missing( experimentAccession ) ) {
        
        warning( "Accession missing. Cannot validate." )

        return( FALSE )
    }

    if( !grepl( "^E-\\w{4}-\\d+$", experimentAccession ) ) {
        
        warning( 
            paste( 
                "\"", 
                experimentAccession, 
                "\" does not look like an ArrayExpress experiment accession. Please check.", 
                sep="" 
            ) 
        )
        
        return( FALSE )

    } else {

        return( TRUE )
    }
}

