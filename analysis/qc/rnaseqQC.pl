#!/usr/bin/env perl
#

use strict;
use warnings;
use 5.10.0;

use Atlas::AtlasConfig::Reader qw( parseAtlasConfig );
use File::Spec;
use Log::Log4perl;
use IPC::Cmd qw( can_run );

$| = 1;

my $logger_config = q(
	log4perl.rootlogger						= INFO, SCREEN
	log4perl.appender.SCREEN				= Log::Log4perl::Appender::Screen
	log4perl.appender.SCREEN.stderr			= 0
	log4j.PatternLayout.cspec.Q				= sub { return "[QC]" }
	log4perl.appender.SCREEN.layout			= Log::Log4perl::Layout::PatternLayout
	log4perl.appender.SCREEN.layout.ConversionPattern = %-5p - %Q %m%n
);

# Initialise logger.
Log::Log4perl::init( \$logger_config );
my $logger = Log::Log4perl::get_logger;

my $atlasProdDir = $ENV{ "ATLAS_PROD" };
unless( $atlasProdDir ) {
	$logger->logdie( "ATLAS_PROD environment variable is not defined, cannot continue." );
}

# Path to script for checking RNA-seq QC results.
my $getQCresultsScript = File::Spec->catfile(
	$atlasProdDir,
	"sw",
	"atlasinstall_prod",
	"atlasprod",
	"irap",
	"single_lib",
	"db",
	"scripts",
	"findCRAMFiles.sh"
);

# Check this user can run the QC results script.
unless( can_run( $getQCresultsScript ) ) {
	$logger->logdie(
		"Cannot run script to get RNA-seq QC results: $getQCresultsScript"
	);
}

my ( $expAcc ) = @ARGV;

my $usage = "Usage:
	rnaseqQC.pl <experiment accession>
";

unless( $expAcc ) { die $usage; }

unless( $expAcc =~ /^E-\w{4}-\d+$/ ) { die $usage; }

my $atlasXMLfile = "$expAcc-configuration.xml";

unless( -r $atlasXMLfile ) {
	$logger->logdie( "Could not read file $atlasXMLfile" );
}

# Parse experiment config.
$logger->info( "Reading experiment config from $atlasXMLfile ..." );
my $experimentConfig = parseAtlasConfig( $atlasXMLfile );
$logger->info( "Successfully read experiment config." );

# Make sure this is an RNA-seq experiment.
unless( $experimentConfig->get_atlas_experiment_type =~ /rnaseq/ ) {
	$logger->logdie( 
		"$expAcc does not look like an RNA-seq experiment: experiment type is ",
		$experimentConfig->get_atlas_experiment_type
	);
}

use Data::Dumper;
print Dumper( $experimentConfig );


# Get a hash of the run accessions from the experiment config; we don't care
# about QC failures of runs that aren't in the experiment config.
my $configRuns = _get_all_config_runs( $experimentConfig );

$logger->info( "Retrieving QC results via $getQCresultsScript ..." );
# Get the RNA-seq QC results
my $rnaseqQCresults = `$getQCresultsScript $expAcc 2>&1`;

# Check whether RNA-seq QC results script ran successfully.
if( $? ) {
	$logger->logdie( "Errors encountered running script $getQCresultsScript\n$rnaseqQCresults" );
}
else {
	$logger->info( "Successfully retrieved QC results." );
}

$logger->info( "Parsing QC results..." );
# Otherwise, what we have should be the table of results for each run in the
# experiment.
my @resultsRows = split /\n/, $rnaseqQCresults;

# A hash to collect failed runs.
my $failedRuns = {};

foreach my $row ( @resultsRows ) {

	# Skip the header.
	if( $row =~ /^study/ ) { next; }

	my @splitRow = split /\t/, $row;

	my $runAcc = $splitRow[ 1 ];
	my $qcStatus = $splitRow[ 5 ];

	# Skip if this run is not in the experiment config.
	unless( $configRuns->{ $runAcc } ) { next; }

	# Add failed run accessions to the hash.
	unless( $qcStatus eq "completed" ) {
		
		$logger->warn( "$runAcc failed QC with status: \"$qcStatus\"" );

		$failedRuns->{ $runAcc } = 1;
	}
}

$logger->info( "Successfully parsed QC results." );

my $qcFileName = $expAcc . "-findCRAMFiles-report.tsv";

$logger->info( "Writing QC results to file $qcFileName ..." );

open( my $fh, ">", $qcFileName ) or $logger->logdie( "Cannot write to file $qcFileName : $!" );
say $fh $rnaseqQCresults;
close $fh;

$logger->info( "Successfully written QC results." );


# If there were any failed runs, go through the XML and remove them.
if( keys %{ $failedRuns } ) {

	# Remove from config.
	$experimentConfig = _remove_rejected_runs( $experimentConfig, $failedRuns );

	# Re-write config now we've removed the failed runs.
	$logger->info( "Renaming original XML config file to \"$atlasXMLfile.beforeQC\"" );
	
	`mv $atlasXMLfile $atlasXMLfile.beforeQC`;
	
	if( $? ) {
		$logger->logdie( "Could not rename file: $!" );
	}
	
	$logger->info( "Writing new XML config file without assays that failed QC..." );
	$experimentConfig->write_xml( "." );
	$logger->info( "Successfully written new XML config file." );

	# Because the Atlas::AtlasConfig modules write files with ".auto" on the end,
	# rename the new one so that it doesn't.
	$logger->info( "Removing \".auto\" from new XML config filename..." );
	`mv $atlasXMLfile.auto $atlasXMLfile`;
	$logger->info( "Successully finished all QC processing." );

} else {
	$logger->info( "All runs passed QC" );
}



sub _get_all_config_runs {

	my ( $experimentConfig ) = @_;

	my $configRuns = {};

	foreach my $analytics ( @{ $experimentConfig->get_atlas_analytics } ) {

		foreach my $assay ( @{ $analytics->get_assays } ) {

			$configRuns->{ $assay->get_name } = 1;
		}
	}

	return $configRuns;
}


sub _remove_rejected_runs {

	my ( $experimentConfig, $failedRuns ) = @_;

	foreach my $runAcc ( keys %{ $failedRuns } ) {

		$logger->info( "Removing run \"$runAcc\" from XML config..." );

		$experimentConfig->remove_assay( $runAcc );

		$logger->info( "Successfully removed \"$runAcc\"" );
	}

	return $experimentConfig;
}
