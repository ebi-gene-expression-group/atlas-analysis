#!/usr/bin/env perl
#
=pod

=head1 NAME

generateFactorsConfig.pl - create an XML factors config file for a baseline Expression Atlas experiment.

=head1 SYNOPSIS

generateFactorsConfig.pl -e E-MTAB-3819 -c E-MTAB-3819-configuration.xml -n "Tissues - Blueprint"

=head1 DESCRIPTION

This script takes an ArrayExpress experiment accession and an Expresison Atlas
experiment configuration file, and generates an XML factors configuration file
for the Expression Atlas user interface. 
    
This file, being automatically generated, will almost always require some
manual curation before it can be used in production.

The display name should follow the same formatting convention as other baseline
expeirments, with the "main" factor type first in Title Case, a " - ", and then
the number of values and the first/main author or consortium name. Note that we
use the word "Tissues" instead of "organism part". 

For example:

 - Tissues - 49 FANTOM5 project
 - Tissues - 39 Roslin
 - Tissues - 5 Palumbo et al
 - Cell Lines - 675 Genentech
 - Developmental Stages - modENCODE

See the baseline experiments page for more examples. If in doubt, check with
the other curators.

=head1 OPTIONS

=over 2

=item -e, --experiment

Required. ArrayExpress experiment accession.

=item -c, --config

Required. Expression Atlas experiment configuration file in XML format.

=item -n, --name

Required. Display name for experiment, to show on
http://www.ebi.ac.uk/gxa/baseline/experiments . See above for examples.

=item -q, --query_factor

Optional. Specify which factor to use as the default query factor.

=item -o, --output_dir

Optional. Where to write the factors config file. If not provided, the current
working directory is used.

=item -d, --debug

Optional. Print additional logging messages, for debugging.

=item -h, --help

Optional. Print this help text.

=back

=head1 AUTHOR

Expression Atlas team <arrayexpress-atlas@ebi.ac.uk>

=cut

use strict;
use warnings;
use 5.10.0;

use Getopt::Long;
use Pod::Usage;
use Log::Log4perl qw( :easy );
use Cwd qw();
use Data::Dumper;
use File::Spec;

use Atlas::Common qw( 
    make_ae_idf_path 
    create_magetab4atlas    
);
use Atlas::AtlasConfig::FactorsConfigFactory qw( create_factors_config );
use Atlas::AtlasConfig::Reader qw( parseAtlasConfig );

my $args = parse_args();

Log::Log4perl->easy_init(
    {
        level   => $args->{ "debug" } ? $DEBUG : $INFO,
        layout  => '%-5p - %m%n',
        file    => "STDOUT",
    }
);

my $logger = Log::Log4perl::get_logger;

# Parse MAGE-TAB and get assays.
my $magetab = create_magetab4atlas( $args->{ "experiment_accession" } );

# Get the assays into a hash.
my %atlasAssays = map { $_->get_name => $_ } @{ $magetab->get_assays };

# Filter out assays that are not found in XML config.
my $atlasAssays = filter_out_nonconfig_assays( \%atlasAssays, $args->{ "experiment_config" } );

# Create new FactorsConfig objects.
my $factorsConfig = create_factors_config( 
    $atlasAssays, 
    $args
);

my $outputFilename = File::Spec->catfile(
    $args->{ "output_directory" },
    $args->{ "experiment_accession" } . "-factors.xml"
);

$factorsConfig->write_xml( $outputFilename );



sub parse_args {

    my %args;

    my $want_help;

    GetOptions(
        "h|help"            => \$want_help,
        "e|experiment=s"    => \$args{ "experiment_accession" },
        "c|config=s"        => \$args{ "experiment_config" },
        "n|name=s"          => \$args{ "display_name" },
        "o|output_dir=s"    => \$args{ "output_directory" },
        "q|query_factor=s"  => \$args{ "default_query_factor" },
        "u|url=s"           => \$args{ "provider_url" },
        "p|provider=s"      => \$args{ "provider_description" },
        "f|fort_lauderdale" => \$args{ "fort_lauderdale" },
        "d|debug"           => \$args{ "debug" }
    );

    unless( $args{ "experiment_accession" } && $args{ "experiment_config" } && $args{ "display_name" } ) {
        
        my $message = 
            "You must specify and experiment accession AND an experiment XML config file name AND a display name.\n";
        
        pod2usage(
            -message    => $message,
            -exitval    => 255,
            -output     => \*STDOUT,
            -verbose    => 1
        );
    }

    unless( $args{ "output_directory" } ) {
        
        say "WARN  - No output directory specifed, will write factors configuration in ", Cwd::cwd();
        $args{ "output_directory" } = Cwd::cwd();
    }

    unless( -w $args{ "output_directory" } ) {
        pod2usage(
            -message => $args{ "output_directory" }. " is not writable or does not exist.\n",
            -exitval => 255,
            -output => \*STDOUT,
            -verbose => 1,
        );
    }

    # Check formatting of displa
    unless( $args{ "display_name" } =~ /^[A-Z].* - .*$/ ) {
        
        say "WARN  - Display name \"" . $args{ "display_name" } . "\" may not be formatted as required.";
        say "WARN  - Please ensure display name follows convention, e.g.:";
        say "WARN  - \"Tissues - 10 Mayer\" or \"Developmental Stages - modENCODE\"";
        say "WARN  - See http://www.ebi.ac.uk/gxa/baseline/experiments for more examples.";
    }

    return \%args;
}


sub filter_out_nonconfig_assays {

    my ( $atlasAssays, $expConfigFile ) = @_;
    
    $logger->debug( "Removing assays not found in XML config..." );

    my $experimentConfig = parseAtlasConfig( $expConfigFile );

    # Check this is a baseline experiment.
    my $experimentType = $experimentConfig->get_atlas_experiment_type;
    
    # Die if not.
    unless( $experimentType =~ /baseline/ ) {
        $logger->logdie( "This is not a baseline experiment. Cannot continue." );
    }

    # Collect the config assay names into a hash.
    my $configAssayNames = {};

    foreach my $analytics ( @{ $experimentConfig->get_atlas_analytics } ) {

        foreach my $assay ( @{ $analytics->get_assays } ) {

            $configAssayNames->{ $assay->get_name } = 1;
        }
    }

    # Go through the magetab assays. If an assay is not found in the experiment
    # config, remove it.
    foreach my $assayName ( keys %{ $atlasAssays } ) {

        unless( $configAssayNames->{ $assayName } ) {

            $logger->debug( 
                $assayName,
                " not found in XML config, removing it."
            );

            delete $atlasAssays->{ $assayName };
        }
    }
    
    $logger->debug( "Finished removing non-config assays." );

    return $atlasAssays;
}
