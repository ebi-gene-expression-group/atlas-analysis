#!/usr/bin/env Rscript

library( methods )

args <- commandArgs( TRUE )

if( length( args ) != 3 ) {
	
	stop( "Usage:\n\n\tcreateAtlasExperimentSummary.R <experiment accession> <atlas experiment directory>" )
}
else {
	source( "/ebi/microarray/home/mkeays/Atlas/git/atlasprod/analysis/atlasExperimentInR/functions.R" )
}

