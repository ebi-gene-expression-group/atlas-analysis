#!/usr/bin/env perl
#

use strict;
use warnings;
use 5.10.0;

use Atlas::AtlasConfig::Reader qw( parseAtlasConfig );
use Atlas::Common qw( create_atlas_site_config );
use File::Spec;
use Log::Log4perl;
use IPC::Cmd qw( can_run );
use Config::YAML;
use Bio::EBI::RNAseqAPI;
use Data::Dumper;

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
	rnaseqQC.pl <experiment accession> <atlas processing directory>
";

unless( $expAcc ) { die $usage; }

unless( $expAcc =~ /^E-\w{4}-\d+$/ ) { die $usage; }

my $atlasProdDir = $ENV{ "ATLAS_PROD" };
unless( $atlasProdDir ) {
	$logger->logdie( "ATLAS_PROD environment variable is not defined, cannot continue." );
}

my $atlasSiteConfig = create_atlas_site_config;

# Path to file listing non-standard experiments. This file contains a list of
# experiments for which there is no QC information stored in the database, and
# hence for which the QC process will not work.
my $nonStandardExpsFile = File::Spec->catfile(
    $atlasProdDir,
    $atlasSiteConfig->get_non_standard_experiments_file
);

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
    
    $logger->warn( "$expAcc is on the list of experiments with no QC information. Skipping QC checks." );

    exit;
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

$logger->info( "Retrieving QC results via RNA-seq API ..." );

# A flag to set if we need to add a new accession to the "too short reads"
# file.
my $newTooShort = 0;

# Collect all the run accessions found in the QC results so that we can check
# that none are missing afterwards. Also remember failed runs, and runs with
# lower than 70% mapped reads.
$_ = {} for my ( $passedQC, $failedRuns, $lowMappedReads, $missing, $inProgress );

# A new RNA-seq API object to run queries with.
my $rnaseqAPI = Bio::EBI::RNAseqAPI->new;

# Put the run accessions from the experiment config into an array for easy
# access.
my @configRunAccs = keys %{ $configRuns };

# Get the results for these runs from the RNA-seq API. Set minimum_mapped_reads
# to zero to get all of the runs.
my $runInfo = $rnaseqAPI->get_runs_by_list( runs => \@configRunAccs, minimum_mapped_reads => 0 );

# Create a hash mapping the run accessions to their info from the API. This
# makes it easier to access their results.
my $mappedRunInfo = _map_info_by_run_acc( $runInfo );

# Write the API we got results into a file for the record. These results can
# change so it's good to record what we got this time.
_write_api_results_to_file( $runInfo, $expAcc );

# Go through the run accessions from the experiment config, and check their status.
foreach my $runAcc ( @configRunAccs ) {

    # Check the run exists in the API results.
    unless( $mappedRunInfo->{ $runAcc } ) {

        $logger->warn( "$runAcc is missing from API results." );

        $missing->{ $runAcc } = 1;

        next;
    }

    # Check if the run is still in progress in ISL.
    if( $mappedRunInfo->{ $runAcc }->{ "STATUS" } eq "In_progress" ) {

        $logger->warn( "$runAcc is still in progress in ISL." );

        $inProgress->{ $runAcc } = 1;

        next;
    }

    # Any status other than Complete is considered a failure from here.
    unless( $mappedRunInfo->{ $runAcc }->{ "STATUS" } eq "Complete" ) {

        $logger->warn( 
            "QC_FAIL: ",
            $runAcc,
            " has failed QC with status \"",
            $mappedRunInfo->{ $runAcc }->{ "STATUS" },
            "\""
        );

        $failedRuns->{ $runAcc } = 1;

        # If the status is about reads being too short, save this.
        if( $mappedRunInfo->{ $runAcc }->{ "STATUS" } =~ /FastqInfo: Read size smaller than/i ) {

            unless( $tooShortAccs->{ $expAcc } ) {
                $newTooShort++;
                $tooShortAccs->{ $expAcc } = 1;
            }
        }

        next;
    }
    
    # If we got here, the run must have passed QC and been processed OK.
    else {

        $passedQC->{ $runAcc } = 1;

        # Store for logging later if percentage mapped reads < 70%.
        if( $mappedRunInfo->{ $runAcc }->{ "MAPPING_QUALITY" } < 70 ) {

            $lowMappedReads->{ $runAcc } = $mappedRunInfo->{ $runAcc }->{ "MAPPING_QUALITY" };
        }
    }
}

# Do some sanity checking.
# We don't want to continue if not all of the runs were present in the API results
if( keys %{ $missing } ) {
    
    my $missingRuns = join "\n", keys %{ $missing };

    $logger->logdie(
        "The following runs are missing from the RNA-seq API results:\n",
        $missingRuns,
        "\nCannot continue when runs are missing from API results."
    );
}

# We don't want to continue if any of the runs are still doing ISL processing.
if( keys %{ $inProgress } ) {

    my $inprogRuns = join "\n", keys %{ $inProgress };

    $logger->logdie(
        "The following runs are still in progress in iRAP single-lib:\n",
        $inprogRuns,
        "\nCannot continue when some runs are still in progress in iRAP single-lib."
    );
}

$logger->info( "All runs in XML config have been QC'ed." );

$logger->info( "Checking mapping quality for runs that passed QC..." );

if( keys %{ $lowMappedReads } ) {

    foreach my $runAcc ( sort keys %{ $lowMappedReads } ) {
        
        $logger->warn( "LOW_MQ: Run ", $runAcc, " has less than 70% reads mapped (", $lowMappedReads->{ $runAcc }, "%)" );
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
    
    my $totalRuns = keys %{ $configRuns };
    my $totalFailed = keys %{ $failedRuns };
    my $totalPassed = $totalRuns - $totalFailed;

    my $pctPassed = ( $totalPassed / $totalRuns ) * 100;

    $logger->info( "PCT_PASSED: $pctPassed % ($totalPassed/$totalRuns) of the runs passed QC." );

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

sub _map_info_by_run_acc {

    my ( $runInfo ) = @_;

    my %mappedRunInfo = map { $_->{ "RUN_IDS" } => $_ } @{ $runInfo };

    return \%mappedRunInfo;
}

sub _write_api_results_to_file {

    my ( $runInfo, $expAcc ) = @_;

    my $apiResultsFile = $expAcc . "-rnaseq-api-results.tsv";

    $logger->info( 
        "Writing API results to ",
        $apiResultsFile,
        " ..."
    );

    open my $fh, ">", $apiResultsFile 
        or $logger->logdie( 
        "Cannot open $apiResultsFile for writing: $!"
    );

    say $fh Dumper( $runInfo );

    close $fh;

    $logger->info(
        "Successfully written API results."
    );
}

