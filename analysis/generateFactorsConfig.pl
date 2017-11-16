#!/usr/bin/env perl
#
=pod

=head1 NAME

generateFactorsConfig.pl - create an XML factors config file for a baseline Expression Atlas experiment.

=head1 SYNOPSIS

generateFactorsConfig.pl -e E-GEOD-26284 -c E-GEOD-26284-configuration.xml --import "annotare" -n "Cell Lines - ENCODE"

generateFactorsConfig.pl -e E-MTAB-2706 -c E-MTAB-2706-configuration.xml -n "Cell Lines - 675 Genentech" -u http://www.gene.com/ -p Genentech

=head1 DESCRIPTION

This script takes an ArrayExpress experiment accession and an Expresison Atlas
experiment configuration file, and generates an XML factors configuration file
for the Expression Atlas user interface. 
    
This file, being automatically generated, will almost always require some
manual curation before it can be used in production.

The display name should follow the same formatting convention as other baseline
expeirments, with the "main" factor type first in Title Case, a " - ", and then
the first/main author or consortium name. Note that we use the word "Tissues"
instead of "organism part". 

For example:

 - Tissues - FANTOM5 project
 - Tissues - Roslin
 - Tissues - Palumbo et al
 - Cell Lines - Genentech
 - Developmental Stages - modENCODE

See the baseline experiments page for more examples. If in doubt, check with
the other curators.

NB: speciesMapping field population is not implemented (yet). This field is
used when the reference genome species is different from the SDRF (RNA sample)
species. For example, data from Pongo pygmaeus may be mapped against the Pongo
abelii genome, because the P. pygmaeus genome is not available and P. abelii is
a closely related species.

Here is an example of how to manually complete the speciesMapping field, if you
need it. The "samples" element is for the SDRF species (note the capial "P"),
while the "genes" element is for the reference genome species. The "genes"
field is used by the UI code to retrieve the correct gene IDs.

 <speciesMapping>
    <genes>pongo abelii</genes>
    <samples>Pongo pygmaeus</samples>
 </speciesMapping>


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

=item -s, --sequence

Optional. Add the flag that the heatmap columns should follow the same sequence
as the assay groups in the experiment XML config file.

=item -u, --url

Optional. Specify the data provider URL, e.g. consortium homepage. This is only used for datasets from
well-known consortia e.g. BLUEPRINT, Genentech. If this option is passed, the
provider name must also be passed, using the "-p" option. The URL and name are
used to display a link to the data provider URL beneath the experiment title on
the Atlas experiment page.

=item -p, --provider

Optional. Specify the data provider name, "The BLUEPRINT Epigenome project". If
this option is passed, the provider URL must also be passed, using the "-u"
option. The URL and name are used to display a link to the data provder beneath
the title on the Atlas experiment page.

=item -a, --agreement

Optional. Specify a data usage agreement be displayed when users download the
data from the Atlas website. E.g. the BLUEPRINT datasets require this.
Currently allowed terms for this option are:

fortLauderdale : BLUEPRINT agreement (E-MTAB-3819, E-MTAB-3827, E-MTAB-4754)

zebrafish : Sanger zebrafish dataset (E-ERAD-475)

=item -o, --output_dir

Optional. Directory to write the factors config file. If not provided, the
current working directory is used.

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
    create_non_strict_magetab4atlas    
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
my $magetab = create_non_strict_magetab4atlas( $args->{ "experiment_accession" }, $args->{ "import_path" } );

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
    $args->{ "experiment_accession" } . "-factors.xml.auto"
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
        "i|import=s"        => \$args{ "import_path" }, # geo, annotare or ena
        "o|output_dir=s"    => \$args{ "output_directory" },
        "q|query_factor=s"  => \$args{ "default_query_factor" },
        "u|url=s"           => \$args{ "provider_url" },
        "p|provider=s"      => \$args{ "provider_name" },
        "a|agreement=s"       => \$args{ "agreement" },
        "s|sequence"        => \$args{ "sequence" },
        "d|debug"           => \$args{ "debug" }
    );

    unless( $args{ "experiment_accession" } && $args{ "experiment_config" } && $args{ "display_name" } && $args{ "import_path" } ) {
        
        my $message = 
            "You must specify and experiment accession AND an experiment XML config file name AND import path AND a display name.\n";
        
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

    if( $args{ "provider_url" } ) {
        
        unless( $args{ "provider_name" } ) {

            pod2usage(
                -message    => "ERROR - Detected data provider URL but no data provider name. Please also pass a data provider name.\n",
                -exitval    => 255,
                -output     => \*STDOUT,
                -verbose    => 1
            );
        }
    }

    if( $args{ "provider_name" } ) {

        unless( $args{ "provider_url" } ) {

            pod2usage(
                -message    => "ERROR - Detected data provider name but no data provider URL. Please also pass a data provider URL.\n",
                -exitval    => 255,
                -output     => \*STDOUT,
                -verbose    => 1
            );
        }
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
