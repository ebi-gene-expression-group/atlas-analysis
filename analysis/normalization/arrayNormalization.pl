#!/usr/bin/env perl
#
# This script parsed MAGE-TAB for a given experiment accession and runs
# microarray normalization by calling an R script.
#
# NB: this script is very similar to arrayQC.pl so maybe (a) they should be
# combined, or (b) we should pull out common parts to a module they both use.

use strict;
use warnings;
use 5.10.0;

# XML config parsing.
use Atlas::AtlasConfig::Reader qw( parseAtlasConfig );
# Site config
use Atlas::Common qw( get_atlas_site_config );
# MAGE-TAB parsing.
use Atlas::Magetab4Atlas;

use File::Spec;
use EBI::FGPT::Config qw( $CONFIG );
use IPC::Cmd qw( can_run );

my $atlasProdDir = $ENV{ "ATLAS_PROD" };
my $atlasSiteConfig = get_atlas_site_config;

# Experiment accession.
my $exptAccession = shift;

unless( $exptAccession ) {
    die "Please provide experiment accession as an argument.\n";
}

# Filename of R script for normalization.
my $normalizationRscript = "arrayNormalization.R";
unless( can_run( $normalizationRscript ) ) {
    die "Script \"$normalizationRscript\" not found. Please ensure it is in your \$PATH and you can run it.\n";
}

# Check that we can run R.
unless( can_run( "R" ) ) {
    die "R was not found. Please ensure it is installed and you can run it.\n";
}

# Path to directory with ArrayExpress/Atlas load directories underneath.
my $exptsLoadStem = File::Spec->catdir( $CONFIG->get_AE2_LOAD_DIR, "EXPERIMENT" );

# miRBase mapped array designs -- we need to subset probes if we find one of these.
# Get an array of miRBase mapping files.
my $miRBaseDirectory = File::Spec->catdir( $atlasProdDir, $atlasSiteConfig->get_mirbase_mappings_directory );
my @A_miRBaseFiles = glob( "$miRBaseDirectory/*.A-*.tsv" );

# Create a hash for easy checking.
my $H_miRBaseFileHash = {};
foreach my $miRBaseFile (@A_miRBaseFiles) {
	# Get the array design from the file name.
	(my $arrayDesign = $miRBaseFile) =~ s/.*(A-\w{4}-\d+)\.tsv/$1/;

	# Add the miRBase mapping file to the hash with the array design as key.
	$H_miRBaseFileHash->{ $arrayDesign } = $miRBaseFile;
}


# Atlas XML config file name.
my $atlasXMLconfigFile = $exptAccession . "-configuration.xml";
# Atlas experiment directory.
my $atlasExperimentDir = File::Spec->catdir( $atlasProdDir, "analysis", "differential", "microarray", "experiments", $exptAccession );
# Full path to XML config file.
my $atlasXMLconfigPath = File::Spec->catfile( $atlasExperimentDir, $atlasXMLconfigFile );

# Parse the config.
my $experimentConfig = parseAtlasConfig( $atlasXMLconfigPath );

# Check that the experiment type is a microarray one.
unless( $experimentConfig->get_atlas_experiment_type =~ /array/ ) {
	die "This does not look like a microarray experiment. Experiment type is \"", $experimentConfig->get_atlas_experiment_type, "\n";
}

# Get the pipeline (e.g. MEXP, MTAB, GEOD, ...) for this experiment.
(my $pipeline = $exptAccession) =~ s/E-(\w{4})-\d+/$1/;

# Experiment load directory and IDF filename.
my $loadDir = File::Spec->catdir( $exptsLoadStem, $pipeline, $exptAccession );
my $idfFilename = File::Spec->catfile( $loadDir, "$exptAccession.idf.txt" );

# Read the MAGE-TAB.
print "Reading MAGE-TAB...\n";
my $magetab4atlas = Atlas::Magetab4Atlas->new( "idf_filename" => $idfFilename );
print "Read MAGE-TAB.\n";

use Data::Dumper;
print Dumper( $magetab4atlas );

# Create hash mapping assay names to raw data file names for each array design:
#
# 	$arraysToAssaysToFiles->{ <array design 1> }->{ <assay 1> } = <file 1>
# 												->{ <assay 2> } = <file 2>
#						  ->{ <array design 2> }->{ <assay 3> } = <file 3>
# 		  				  ...
# Also pass $experimentType and get back normalization mode to pass to R script.
my ($H_arraysToAssaysToFiles, $normalizationMode) = &makeArraysToAssaysToFiles( $magetab4atlas, $loadDir, $experimentConfig);

# Log how many array designs and assays were found
my $arrayCount = keys %{ $H_arraysToAssaysToFiles };
print "Found $arrayCount array design";
if($arrayCount > 1) { print "s"; }
print ":\n";
foreach my $arrayDesign (keys %{ $H_arraysToAssaysToFiles }) {
	my $assayCount = keys %{ $H_arraysToAssaysToFiles->{ $arrayDesign }};
	print "\t$arrayDesign: $assayCount assays.\n";
}
# Also what kind of normalization will be done. This is either:
# 	- "oligo" : using oligo package for Affymetrix arrays.
# 	- "agil1" : using limma pacakge for Agilent 1-colour data.
# 	- "agil2" : using limma package for Agilent 2-colour data.
print "The normalization mode is $normalizationMode.\n";

# Run normalization for each array design.
foreach my $arrayDesign (keys %{ $H_arraysToAssaysToFiles }) {
	# Check if the array design has a miRBase mapping file
	# Flag
	my $miRBaseFile = 0;
	if(exists($H_miRBaseFileHash->{ $arrayDesign })) {
		print "$arrayDesign is a microRNA array design.\n";
		$miRBaseFile = $H_miRBaseFileHash->{ $arrayDesign };
	}

	# Write a file to read into R
	my $tempFile = "/tmp/$exptAccession"."_$arrayDesign.$$.tsv";
	open(my $tmpFH, ">", $tempFile) or die "Can't create file \"$tempFile\": $!\n";

	# Write headers
	print $tmpFH "AssayName\tFilename";
	foreach my $assayName (keys %{ $H_arraysToAssaysToFiles->{ $arrayDesign }}) {
		print $tmpFH "\n$assayName\t$H_arraysToAssaysToFiles->{ $arrayDesign }->{ $assayName }";
	}
	close $tmpFH;

	# Name for normalized data file.
	my $normalizedDataFile = $exptAccession."_".$arrayDesign."-normalized-expressions.tsv.undecorated";

	print "Running normalization in R for $exptAccession, array design $arrayDesign...\n";

	# Run R script to do normalization with Bioconductor packages in R.
	# NB: Using 2>&1 means that nothing from R is printed to STDOUT. If there's an
	# error, then it's printed by the part following this line.
	my $RscriptOutput = `$normalizationRscript $tempFile $normalizationMode $normalizedDataFile $miRBaseFile 2>&1`;
	
	# If there was an error in R die and print the error.
	if($RscriptOutput =~ /error/i) {
		die "Error encountered during normalization of $exptAccession on array $arrayDesign. Full output from R is below:
		------------------------
		$RscriptOutput
		";
	}
	else {
		# For 2-colour data, rename files created.
		if($normalizationMode eq "agil2") {
			my $logFCfile = $exptAccession."_$arrayDesign-log-fold-changes.tsv.undecorated";
			my $aValuesFile = $normalizedDataFile.".A-values";
			my $avgIntensitiesFile = $exptAccession."_$arrayDesign-average-intensities.tsv.undecorated";

			`mv $normalizedDataFile $logFCfile`;
			`mv $aValuesFile $avgIntensitiesFile`;
		}

		print "Normalization for array design $arrayDesign completed.\n";
	}

	# Delete temporary file.
	`rm $tempFile`;
}


# Subroutines

# &makeArraysToAssaysToFiles
# 
# 	- Creates a hash matching each data file to each assay, for each array design.
# 	- E.g.:
# 		$H->{ <array design 1> }->{ <assay 1> } = <file 1>
# 								->{ <assay 2> } = <file 2>
# 		  ->{ <array design 2> }->{ <assay 3> } = <file 3>
# 		  ...
# Arguments:
# 	- $magetab4atlas : a Atlas::Magetab4Atlas object.
# 	- $loadDir : path to load directory containing raw data files.
# 	- $experimentConfig : Atlas::AtlasConfig::ExperimentConfig object.
sub makeArraysToAssaysToFiles {
	# Atlas::Magetab4Atlas object and path to load directory.
	my ($magetab4atlas, $loadDir, $experimentConfig ) = @_;
	
	my $experimentType = $experimentConfig->get_atlas_experiment_type;

	# Get all the assay names from the experiment config.
	my $allAnalytics = $experimentConfig->get_atlas_analytics;

	# Create a hash of all the assay names found in the experiment config,
	# mapped to their respective array designs.
	my $arrayDesignsToAssayNames = {};

	foreach my $analytics ( @{ $allAnalytics }) {

		my $analyticsAssays = $analytics->get_assays;

		my $arrayDesign = $analytics->get_platform;

		$arrayDesignsToAssayNames->{ $arrayDesign } = [];

		foreach my $assay ( @{ $analyticsAssays } ) {

			push @{ $arrayDesignsToAssayNames->{ $arrayDesign } }, $assay->get_name;
		}
	}

	# Ref to empty hash to fill
	my $H_arraysToAssaysToFiles = {};
	
	# Normalization mode.
	my $normalizationMode = 0;

	# Go through the assays...
	foreach my $assay4atlas (@{ $magetab4atlas->get_assays }) {

		# Get assay name
		my $assayName = $assay4atlas->get_name;

		# Escape any metacharacters
		my $assayNameEsc = quotemeta( $assayName );
		
		# Array design
		my $arrayDesign = $assay4atlas->get_array_design;

		# Check that we saw this assay in the XML config. If not, skip it.
		unless( grep { /^$assayNameEsc$/ } @{ $arrayDesignsToAssayNames->{ $arrayDesign } } ) { 
			print "Assay \"$assayName\" not found in XML config, not including in normalization.\n";
			next;
		}

		# Raw data filename.
		my $arrayDataFile = File::Spec->catfile( $loadDir, $assay4atlas->get_array_data_file );
		
		# For 1-colour array data, need to tell the R script whether this is
		# Affymetrix or Agilent (or other -- for now we only handle Affy and Agil
		# data). Work this out based on the file extension for now -- another
		# possibility might be to parse the ADF, but there is no standard way
		# to record manufacturer there so that might not be ideal.
		# For now only allow one experiment type per experiment but may
		# need to change this.
		if( $experimentType =~ /1colour/ ) {
			if( $arrayDataFile =~ /\.cel$/i ) { $normalizationMode = "oligo"; }
			else { $normalizationMode = "agil1"; }
		} elsif( $experimentType =~ /2colour/ ) {
			$normalizationMode = "agil2";

			# Remove label name from assay name which was added by Atlas::Magetab4Atlas.
			$assayName =~ s/\.Cy\d$//;
		}
		
		# Add data to hash.
		$H_arraysToAssaysToFiles->{ $arrayDesign }->{ $assayName } = $arrayDataFile;
	}
	
	return ($H_arraysToAssaysToFiles, $normalizationMode);
}
