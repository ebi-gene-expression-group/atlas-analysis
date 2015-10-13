#!/usr/bin/env perl
#

use strict;
use warnings;
use 5.10.0;

use JSON::Parse qw( parse_json );
use Data::Dumper;
use File::Spec;
use IPC::Cmd qw( can_run );
use Log::Log4perl;
use Getopt::Long;
use Cwd qw();

my $logger_config = q(
    log4perl.rootlogger                  = INFO, SCREEN
    log4perl.appender.SCREEN             = Log::Log4perl::Appender::Screen
    log4perl.appender.SCREEN.stderr      = 0
    log4perl.appender.SCREEN.layout      = Log::Log4perl::Layout::PatternLayout
    log4perl.appender.SCREEN.layout.ConversionPattern = %-5p - %m%n
);

Log::Log4perl::init( \$logger_config );
my $logger = Log::Log4perl::get_logger;

my $percentileCalcScript = "calculate_percentile_ranks.R";

unless( can_run( $percentileCalcScript ) ) {
    $logger->logdie( "Cannot run R script to do percentile ranks calculation. Please check that it is in your path." );
}

# Get the experiment accession from the XML config filename.
my $configFilename = ( glob( "*-configuration.xml" ) )[ 0 ];
unless( $configFilename ) {
    $logger->logdie( "No XML configuration file found in ", Cwd::cwd(), " . Cannot determine experiment accession." );
}

# Get the experiment accession from the config filename.
( my $expAcc = $configFilename ) =~ s/^(E-\w{4}-\d+).*/$1/;

# Delete previous files (if any) as they confuse things.
my @previousPercentileRanksFiles = glob( "*-percentile-ranks.tsv" );

if( @previousPercentileRanksFiles ) {
    
    $logger->info( "Removing previous percentile ranks files..." );
    
    `rm *-percentile-ranks.tsv`;

    $logger->info( "Done" );
}

# File the analytics file(s).
my @analyticsFiles = glob( "*-analytics.tsv.unrounded" );

# Check we got at least one analytics filename.
if( @analyticsFiles == 0 ) {
    $logger->logdie( "No unrounded analytics files found in ", Cwd::cwd() );
}

# Go through the analytics files and calculate percentile ranks for each of
# them.
foreach my $analyticsFile ( @analyticsFiles ) {
    if( -r $analyticsFile ) {
        
        $logger->info( "Calculating percentile ranks for $analyticsFile ..." );

        my $r_output = `$percentileCalcScript $analyticsFile 2>&1`;

        if( $r_output =~ /error/i ) {
            $logger->logdie( "Problems during R script run, output from R below.\n--------\n$r_output\n" );
        }

        $logger->info( "Percentile ranks calculation successful." );
    }
    else {
        $logger->logdie( "Cannot read analytics for $expAcc." );
    }
}
    

# Change the filename of the percentile ranks file so that it only has the
# experiment accession. Combine files from multi-array experiments.
my @percentileRankFiles = glob( "*percentile-ranks.tsv" );

# New filename.
my $allRanksFile = $expAcc . "-percentile-ranks.tsv";

# If we have more than one array design in a microarray experiment, there will
# be more than one percentile ranks file. Combine them.
if( @percentileRankFiles > 1 ) {

    $logger->info( "Found multiple rank files for this experiment, combining them..." );

    my $gene2contrast2rank = {};

    my @allContrastIDs = ();

    foreach my $rankFile ( @percentileRankFiles ) {

        open my $fh, "<", $rankFile or die( "Cannot open rank file: $!\n" );
        
        my @contrastIDs;

        while( defined( my $line = <$fh> ) ) {
            
            chomp $line;
            
            my @lineSplit = split( "\t", $line );

            if( $line =~ /^Gene.ID/ ) {

                shift @lineSplit;

                @contrastIDs = @lineSplit;

                push @allContrastIDs, @contrastIDs;
            }
            else {

                my $geneID = shift @lineSplit;

                for( my $i = 0; $i < @lineSplit; $i++ ) {
                    
                    my $contrastID = $contrastIDs[ $i ];
                    my $rank = $lineSplit[ $i ];

                    $gene2contrast2rank->{ $geneID }->{ $contrastID } = $rank;
                }
            }
        }
        close $fh;
    }

    # Write out a new file with the combined contrasts in.

    open( my $fh, ">", $allRanksFile ) or die( "Cannot open $allRanksFile for writing: $!\n" );
    
    # Write header.
    my $joinedContrastIDs = join "\t", @allContrastIDs;
    say $fh "Gene.ID\t$joinedContrastIDs";

    # Write the rest.
    foreach my $geneID ( sort keys %{ $gene2contrast2rank } ) {

        my @ranks2write = ();

        foreach my $contrastID ( @allContrastIDs ) {
            
            my $rank = $gene2contrast2rank->{ $geneID }->{ $contrastID };

            if( defined( $rank ) ) {
                push @ranks2write, $rank;
            } else {
                push @ranks2write, "NA";
            }
        }

        my $joinedRanks = join "\t", @ranks2write;

        say $fh "$geneID\t$joinedRanks";
    }

    close $fh;

    # Deleted the now un-needed original files.
    $logger->info( "Removing files for each array design..." );
    foreach my $originalFile ( @percentileRankFiles ) {
        `rm $originalFile`;
    }
    $logger->info( "Done." );
}
# Otherwise, if this is a microarray exp, just change the filename of the
# original file to remove the array design accession.
elsif( $percentileRankFiles[ 0 ] =~ /A-\w{4}-\d+/ ) {
    `mv $percentileRankFiles[ 0 ] $allRanksFile`;
}

$logger->info( "Successfully calculated percentile ranks for $expAcc." );


