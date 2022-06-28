#!/usr/bin/env perl

use strict;
use warnings;
use 5.10.0;

use Atlas::AtlasConfig::Reader qw( parseAtlasConfig );
use Atlas::Common qw(
    get_supporting_file
    create_atlas_site_config
    check_privacy_in_ae
    get_ena_study_id
);
use File::Spec;
use Log::Log4perl;
use IPC::Cmd qw( can_run );
use Config::YAML;

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

my ( $expAcc, $atlasProcessingDir ) = @ARGV;

my $usage = "Usage:
	rnaseqQC.pl <experiment accession>
";

unless( $expAcc ) { die $usage; }

unless( $expAcc =~ /^E-\w{4}-\d+$/ ) { die $usage; }

my $atlasProdDir = $ENV{ "ATLAS_PROD" };
unless( $atlasProdDir ) {
	$logger->logdie( "ATLAS_PROD environment variable is not defined, cannot continue." );
}

# Path to file listing non-standard experiments. This file contains a list of
# experiments for which there is no QC information stored in the database, and
# hence for which the QC process will not work.
my $nonStandardExpsFile = get_supporting_file("non_standard_experiments.yml");

unless( -f $nonStandardExpsFile ) {
    $logger->logdie(
        "Cannot find file $nonStandardExpsFile"
    );
}

# Check that the accession is not in the list of experiments with no QC info.
my $nonStandardExps = Config::YAML->new(
    config => $nonStandardExpsFile
);

if( grep $_ eq $expAcc, @{ $nonStandardExps->get_no_qc_info } ) {

    $logger->info( "$expAcc is on the list of experiments with no QC information. Skipping QC checks." );

    exit;
}

my $atlasSiteConfig = create_atlas_site_config;

my $irapSingleLib = $ENV{ "IRAP_SINGLE_LIB" };

# Path to script for checking RNA-seq QC results.
my $islResultsScript = File::Spec->catfile(
	$irapSingleLib,
	$atlasSiteConfig->get_isl_results_script
);

# Check this user can run the QC results script.
unless( can_run( $islResultsScript ) ) {
	$logger->logdie(
		"Cannot run script to get RNA-seq QC results: $islResultsScript"
	);
}

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

$logger->info( "Retrieving QC results via $islResultsScript ..." );
# Get the RNA-seq QC results
my $rnaseqQCresults;

$rnaseqQCresults = `$islResultsScript $expAcc 2>&1`;

# Check whether RNA-seq QC results script ran successfully.
if( $? ) {
	$logger->logdie( "Errors encountered running script $islResultsScript\n$rnaseqQCresults" );
}
else {
	$logger->info( "Successfully retrieved QC results." );
}

my $qcFileName = $expAcc . "-irap-single-lib-report.tsv";

$logger->info( "Writing QC results to file $qcFileName ..." );

open( my $qcfh, ">", $qcFileName ) or $logger->logdie( "Cannot write to file $qcFileName : $!" );
say $qcfh $rnaseqQCresults;
close $qcfh;

$logger->info( "Successfully written QC results." );

$logger->info( "Parsing QC results..." );
# Otherwise, what we have should be the table of results for each run in the
# experiment.
my @resultsRows = split /\n/, $rnaseqQCresults;

# Collect all the run accessions found in the QC results so that we can check
# that none are missing afterwards. Also remember failed runs, and runs with
# lower than 70% mapped reads.
$_ = {} for my ( $qcRuns, $passedQC, $failedRuns, $lowMappedReads, $inProgress );

foreach my $row ( @resultsRows ) {

	# Skip the header.
	if( $row =~ /^study/ ) { next; }

	my @splitRow = split /\t/, $row;

	my $runAcc = $splitRow[ 1 ];
	my $qcStatus = $splitRow[ 6 ];

    # Save the run accession.
    $qcRuns->{ $runAcc } = 1;

	# Skip if this run is not in the experiment config.
	unless( $configRuns->{ $runAcc } ) { next; }

	# Add failed run accessions to the hash.
	if( $qcStatus eq "completed" || $qcStatus eq "on_ftp" || $qcStatus eq "completed_no_transcript" || $qcStatus eq 'failed_ftp_copy' || $qcStatus eq 'analysis_complete' ) {

        $passedQC->{ $runAcc } = 1;

        # Also see if we got less than 70% reads mapped. If so, add them to
        # the hash to save them.
        my $percentMappedReads = $splitRow[ 7 ];

        if( $percentMappedReads < 70 ) {

            $lowMappedReads->{ $runAcc } = $percentMappedReads;
        }
    }
    elsif( $qcStatus eq "irap_single_lib" ) {

        $inProgress->{ $runAcc } = 1;
    }
    else {

        $logger->info( "QC_FAIL: $runAcc failed QC with status: \"$qcStatus\"" );

		$failedRuns->{ $runAcc } = 1;
	}
}

$logger->info( "Successfully parsed QC results." );

$logger->info( "Checking that all runs in experiment config have been QC'ed..." );
# Check that all runs from the XML config were found in the QC results. Die if
# not.
my $runsMissing = 0;

foreach my $runAcc ( keys %{ $configRuns } ) {

    unless( $qcRuns->{ $runAcc } ) {

        $logger->info( "QC_FAIL: $runAcc was not found in QC results." );

        $runsMissing++;
    }
}

if( $runsMissing ) {
    $logger->logdie( "Not all runs in XML config were found in QC results. Cannot continue." );
}

$logger->info( "All runs in XML config have been QC'ed." );

# Check if any runs were in progress, die here if so, as the experiment should
# not have been submitted to ISL.
if( keys %{ $inProgress } ) {

    my $inprogRuns = join "\n", keys %{ $inProgress };

    $logger->logdie(
        "The following runs are still in progress in iRAP single-lib:\n",
        $inprogRuns,
        "\nCannot continue when some runs are still in progress in iRAP single-lib."
    );
}

$logger->info( "Checking mapping quality for runs that passed QC..." );

if( keys %{ $lowMappedReads } ) {

    foreach my $runAcc ( sort keys %{ $lowMappedReads } ) {

        $logger->info( "LOW_MQ: Run ", $runAcc, " has less than 70% reads mapped (", $lowMappedReads->{ $runAcc }, "%)" );
    }
}
else {
    $logger->info( "All runs have >= 70% reads mapped." );
}

# Also check that all QC'ed runs exist in the raw counts matrix.
my $countsMatrixFile = File::Spec->catfile( $atlasProcessingDir, $expAcc . "-raw-counts.tsv.undecorated" );

$logger->info( "Checking that all runs that passed QC exist in counts matrix $countsMatrixFile" );

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

    my $missingRunString = join "\n", ( keys %{ $missingFromCounts } );

    $logger->logdie( "The following run(s) were not found in the counts matrix:\n$missingRunString" );
}

# If we're still alive, log that all was OK.
$logger->info( "All runs that passed QC were found in the counts matrix." );

# Report what percentage of runs passed QC.
my $totalRuns = keys %{ $configRuns };
my $totalFailed = keys %{ $failedRuns };
my $totalPassed = $totalRuns - $totalFailed;

my $pctPassed = ( $totalPassed / $totalRuns ) * 100;

$logger->info( "PCT_PASSED: $pctPassed % ($totalPassed/$totalRuns) of the runs passed QC." );

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