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

    # Only use valid accessions to download.
    experimentAccessions <- experimentAccessions[ 
        which( 
            sapply( 
                experimentAccessions, function( accession ) {
                    .isValidExperimentAccession( accession )
                }
            )
        )
    ]

    # The experimentAccessions vector is empty if none of the accessions are
    # valid. Just quit here if so.
    if( length( experimentAccessions ) == 0 ) {
        stop( "None of the accessions passed are valid ArrayExpress accessions. Cannot continue." )
    }
    
    # Go through each one and download it, creating a list.
    # So that the list has the experiment accessions as the names, use them as
    # names for the vector before starting to create the list.
    names( experimentAccessions ) <- experimentAccessions

    experimentSummaryList <- SimpleList( 
        
        lapply( experimentAccessions, function( experimentAccession ) {

            getAtlasExperiment( experimentAccession )
        } 
    ) )

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


# searchAtlasExperiments
#   - Search (currently against ArrayExpress API) for datasets in Atlas matching given terms.
#   - TODO: since we are currently using the ArrayExpress API, we can't query
#   for genes, only sample properties.
searchAtlasExperiments <- function( properties, species = NULL ) {

    if( missing( properties ) ) {
        stop( "Please provide at least one search term." )
    }

    if( typeof( properties ) != "character" ) {
        stop( "Please provide search term(s) as a character vector." )
    }
    
    properties <- sapply( properties, URLencode )

    if( missing( species ) ) {
    
        cat( "No species was provided. Will search for data from all available species.\n" )
    
    } else if( typeof( species ) != "character" ) {
       
       stop( "Please provide species as a character vector." )
    
    } else if( length( species ) > 1 ) {

        stop( "More than one species found. You may only specify one species at a time." )
    }

    aeAPIbase <- "http://www.ebi.ac.uk/arrayexpress/xml/v2/experiments?keywords="

    queryURL <- paste( 
        aeAPIbase,
        paste( properties, collapse = "+OR+" ),
        "&gxa=TRUE",
        sep = ""
    )

    if( !missing( species ) ) {
    
        species <- URLencode( species )

        queryURL <- paste( 
            queryURL, 
            "&species=", 
            species, 
            sep = "" 
        )
    }

    cat( 
        paste( 
            "Searching ArrayExpress for experiments in Atlas using the following URL:\n",
            queryURL,
            " ...\n",
            sep = ""
        )
    )
    
    # Run the query and download the resulting XML.
    aeResultsXmlTree <- xmlInternalTreeParse( queryURL, isURL = TRUE )
    
    # TODO: make sure the above worked i.e. wrap it in try().

    # Pull out the root node ("experiments").
    allExpsNode <- xmlRoot( aeResultsXmlTree )

    # Get a list of all the experiments.
    allExperiments <- xmlElementsByTagName( allExpsNode, "experiment" )

    # TODO: log the number of experiments found.

    # Pull out the title, accession, and species of each experiment.
    # TODO: What else might be useful here?
    #   - experiment type i.e. microarray or RNA-seq
    #   - description(?), 
    #   - factor(s)? 
    #   - ...?
    resultsList <- lapply( allExperiments, function( experimentNode ) {

        expAcc <- xmlValue( xmlElementsByTagName( experimentNode, "accession" )$accession )

        expTitle <- xmlValue( xmlElementsByTagName( experimentNode, "name" )$name )

        species <- xmlValue( xmlElementsByTagName( experimentNode, "organism" )$organism )

        expType <- xmlValue( xmlElementsByTagName( experimentNode, "experimenttype")$experimenttype )

        list( accession = expAcc, title = expTitle, species = species, expType = expType )

    } )
    
    allAccessions <- sapply( resultsList, function( x ) { x$accession } )
    allExpTypes <- sapply( resultsList( function( x ) { x$expType } )
    allSpecies <- sapply( resultsList, function( x ) { x$species } )
    allTitles <- sapply( resultsList, function( x ) { x$title } )

    resultsSummary <- data.frame( 
        Accession = allAccessions, 
        Species = allSpecies, 
        Type = allExpTypes, 
        Title = allTitles, 
        stringsAsFactors = FALSE 
    )
    resultsSummary <- resultsSummary[ order( resultsSummary$Species, resultsSummary$Accession ), ]
    
    return( resultsSummary )
}
