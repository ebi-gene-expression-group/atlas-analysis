#!/usr/bin/env perl
#

use strict;
use warnings;
use 5.10.0;

use Atlas::AtlasConfig::Reader qw( parseAtlasConfig );
use Atlas::Common qw( get_atlas_site_config );
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

my $atlasSiteConfig = get_atlas_site_config;

# Path to script for checking RNA-seq QC results.
my $getQCresultsScript = File::Spec->catfile( 
	$atlasProdDir,
	$atlasSiteConfig->get_find_cram_files_script
);

# Check this user can run the QC results script.
unless( can_run( $getQCresultsScript ) ) {
	$logger->logdie(
		"Cannot run script to get RNA-seq QC results: $getQCresultsScript"
	);
}

# Get the list of accessions of experiments that have runs rejected by iRAP for
# having too-short reads (less than 30bp).
my $tooShortReadsFile = File::Spec->catfile(
	$atlasProdDir,
	$atlasSiteConfig->get_exps_with_runs_rejected_for_small_read_size
);

# Read the accessions of experiments that have been rejected for having
# too-short reads.
my $tooShortAccs = {};
open my $fh, "<", $tooShortReadsFile
	or $logger->logdie( "Cannot read file $tooShortReadsFile: $!" );
while( defined( my $line = <$fh> ) ) {

	chomp $line;
	
	if( $line =~ /^E-\w{4}-\d+$/ ) {
		$tooShortAccs->{ $line } = 1;
	}
}
close $fh;

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

# A flag to set if we need to add a new accession to the "too short reads"
# file.
my $newTooShort = 0;

# Collect all the run accessions found in the QC results so that we can check
# that none are missing afterwards.
$_ = {} for my ( $qcRuns, $passedQC );

foreach my $row ( @resultsRows ) {

	# Skip the header.
	if( $row =~ /^study/ ) { next; }

	my @splitRow = split /\t/, $row;

	my $runAcc = $splitRow[ 1 ];
	my $qcStatus = $splitRow[ 5 ];

    # Save the run accession.
    $qcRuns->{ $runAcc } = 1;

	# Skip if this run is not in the experiment config.
	unless( $configRuns->{ $runAcc } ) { next; }

	# Add failed run accessions to the hash.
	unless( $qcStatus eq "completed" || $qcStatus eq "on_ftp" ) {
		
		$logger->warn( "$runAcc failed QC with status: \"$qcStatus\"" );

		$failedRuns->{ $runAcc } = 1;
		
		if( $qcStatus =~ /FastqInfo: Read size smaller than/i ) {

			unless( $tooShortAccs->{ $expAcc } ) {
				
				$newTooShort++;

				$tooShortAccs->{ $expAcc } = 1;
			}
		}
	}
    # If the run passed QC, save it so we can make sure it exists in the counts
    # matrix.
    else {

        $passedQC->{ $runAcc } = 1;
    }
}

$logger->info( "Successfully parsed QC results." );

$logger->info( "Checking that all runs in experiment config have been QC'ed..." );
# Check that all runs from the XML config were found in the QC results. Die if
# not.
my $runsMissing = 0;

foreach my $runAcc ( keys %{ $configRuns } ) {

    unless( $qcRuns->{ $runAcc } ) {

        $logger->error( "$runAcc was not found in QC results." );
 
        $runsMissing++;
    }
}

if( $runsMissing ) {
    $logger->logdie( "Not all runs in XML config were found in QC results. Cannot continue." );
}

$logger->info( "All runs in XML config have been QC'ed." );

# Also check that all QC'ed runs exist in the raw counts matrix.
$logger->info( "Checking that all runs that passed QC exist in counts matrix..." );

my $countsMatrixFile = File::Spec->catfile( $ENV{ "IRAP_SINGLE_LIB" }, "studies", $expAcc, "genes.raw.tsv" );

my %countsMatrixRuns;

open( my $countsFH, "<", $countsMatrixFile ) or $logger->logdie( "Cannot open $countsMatrixFile: $!" );

while( defined( my $line = <$countsFH> ) ) {

    chomp $line;

    my @headers = split "\t", $line;

    # Remove first element (e.g. "Gene ID") as we don't want it.
    shift @headers;

    # Add remaining header (run accessions) to hash.
    %countsMatrixRuns = map { $_ => 1 } @headers;

    # Quit as we only care about the first line of this file.
    last;
}

close $countsFH;

# Collect missing run accessions to report later.
my $missingFromCounts = {};

# Make sure every run that passed QC exists in the counts matrix.
foreach my $runAcc ( keys %{ $passedQC } ) {

    unless( $countsMatrixRuns{ $runAcc } ) {
        $missingFromCounts->{ $runAcc } = 1;
    }
}

if( keys %{ $missingFromCounts } ) {

    my $missingRunString = join ", ", ( keys %{ $missingFromCounts } );

    $logger->logdie( "The following run(s) were not found in the counts matrix: $missingRunString" );
}

# If we're still alive, log that all was OK.
$logger->info( "All runs that passed QC were found in the counts matrix." );

my $qcFileName = $expAcc . "-findCRAMFiles-report.tsv";

$logger->info( "Writing QC results to file $qcFileName ..." );

open( my $qcfh, ">", $qcFileName ) or $logger->logdie( "Cannot write to file $qcFileName : $!" );
say $qcfh $rnaseqQCresults;
close $qcfh;

$logger->info( "Successfully written QC results." );

# If there were any new accessions to add to the "too-short reads" file,
# re-write it now.
if( $newTooShort ) {

	$logger->info( "Adding experiment accession to $tooShortReadsFile ..." );

	open my $tsfh, ">", $tooShortReadsFile
		or $logger->logdie( "Cannot write to file $tooShortReadsFile : $!" );

	foreach my $acc ( keys %{ $tooShortAccs } ) {
		say $tsfh $acc;
	}

	close( $tsfh );

	$logger->info( "Successfully written accessions." );
}


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

# end.



##############
# Subroutines.

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
