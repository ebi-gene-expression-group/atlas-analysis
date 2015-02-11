# createAtlasSiteConfig
#   - Return a new AtlasSiteConfig object
createAtlasSiteConfig <- function( ) {

    # We need the path of the YAML config file.
    # We can deduce this using the output of the R.home() function, which is
    # the directory where R is installed.
    rHome <- R.home()

    # Split the path on "/" characters so that we can remove everything after
    # the atlasinstall directory more easily. Use the platform's file separator
    # (in UNIX a "/").
    splitPath <- strsplit( rHome, .Platform$file.sep )[[ 1 ]]

    # Find the index of the atlasinstall_<instance> directory.
    atlasInstallIndex <- grep( "atlasinstall", splitPath )

    # Check that we got an index. If not, we must not be using an atlasinstall R.
    if( length( atlasInstallIndex ) == 0 ) {
        stop( "Could not find atlasinstall directory in RHOME path. Unable to deduce path to site config file. Please ensure you are using an atlasinstall R installation." )
    }

    # Remove everything after the atlasinstall directory from the path.
    atlasInstallSplitPath <- head( splitPath, n = atlasInstallIndex )

    # Stick the path back together again.
    atlasInstallPath <- paste( atlasInstallSplitPath, collapse=.Platform$file.sep )

    # Add the path to the site config file.
    atlasSiteConfigFile <- file.path( atlasInstallPath, "atlasprod", "supporting_files", "AtlasSiteConfig.yml" )

    # Load the yaml package.
    suppressMessages( library( yaml ) )

    # Read in the site config file.
    atlasSiteConfig <- yaml.load_file( atlasSiteConfigFile )

    return( atlasSiteConfig )
}
