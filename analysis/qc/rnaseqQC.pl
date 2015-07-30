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

my $atlasXMLfile = "$expAcc-configuration.xml";

unless( -r $atlasXMLfile ) {
	$logger->logdie( "Could not read file $atlasXMLfile" );
}

# Parse experiment config.
my $experimentConfig = parseAtlasConfig( $atlasXMLfile );

# Make sure this is an RNA-seq experiment.
unless( $experimentConfig->get_atlas_experiment_type =~ /rnaseq/ ) {
	$logger->logdie( 
		"$expAcc does not look like an RNA-seq experiment: experiment type is ",
		$experimentConfig->get_atlas_experiment_type
	);
}

# Get the RNA-seq QC results
my $rnaseqQCresults = `$getQCresultsScript $expAcc 2>&1`;

# Check whether RNA-seq QC results script ran successfully.
if( $? ) {
	$logger->logdie( "Errors encountered running script $getQCresultsScript\n$rnaseqQCresults" );
}

# Otherwise, what we have should be the table of results for each run in the
# experiment.
my @resultsRows = split /\n/, $rnaseqQCresults;

foreach my $row ( @resultsRows ) {

	my @splitRow = split /\t/, $row;

	my $runAcc = $splitRow[ 1 ];
	my $qcStatus = $splitRow[ 5 ];

	unless( $qcStatus eq "completed" ) {
		say "$runAcc has status $qcStatus";
	}
}
