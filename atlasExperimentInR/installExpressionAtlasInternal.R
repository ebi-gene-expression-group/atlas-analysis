#!/usr/bin/env Rscript

# Script to check for installation of ExpressionAtlas R package 
# and install if needed

suppressMessages( library( optparse ) )
suppressMessages( library( devtools ) )

# get parsed arguments
args <- parse_args(OptionParser(option_list= list(
	make_option(
		c("-v", "--version"),
		help="version require"
	)
)))

# version of internal Expression Atlas R package in use
require_version<-args$version

# check for the current version of ExpressionAtlasInternal is similar to required version
# if not equal then install the required version
if ( packageDescription("ExpressionAtlasInternal")$Version != require_version ) {
    cat("ExpressionAtlasInternal  version", require_version, "not installed.\n")

    cat("Installing ExpressionAtlasInternal version", require_version,"\n")
    devtools::install_github( repo="ebi-gene-expression-group/internal-ExpressionAtlas", ref = require_version)
    
    # secondary check post installation of the package to validate the version tag what we have installed is also the same in the DESCRIPTION file.
    if (packageDescription("ExpressionAtlasInternal")$Version != require_version){
		stop("Tag version is not the same as in DESCRIPTION file. Make sure DESCRIPTION version ", packageDescription("ExpressionAtlasInternal")$Version, " and tag ", require_version, " are the same\n")
	}

} else {
	  cat("Required ExpressionAtlasInternal version", require_version, "is in use\n")
}
