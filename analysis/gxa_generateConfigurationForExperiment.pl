#!/usr/bin/env perl
#
=pod

=head1 NAME

gxa_generateConfigurationForExperiment.pl - create an XML config file for an Expression Atlas experiment.

=head1 SYNOPSIS

gxa_generateConfigurationForExperiment.pl -e E-MTAB-1066 -t differential

gxa_generateConfigurationForExperiment.pl -e E-MTAB-513 -t baseline -l paired

=head1 DESCRIPTION

This script takes an ArrayExpress experiment accession and an analysis type
(baseline or differential), and creates an XML config file for Expression
Atlas.

=head1 OPTIONS

=over 2

=item -e --experiment

Requred. ArrayExpress accession of experiment.

=item -t --type

Required. Type of analysis. Must be one of "baseline" or "differential".

=item -l --library

Optional. Specify the type of RNA-seq libraries to retrieve from MAGE-TAB. Must
be one of "paired" or "single".

=item -r --reference

Optional. Differential experiments only. Specify a value to use as the
reference factor value in contrasts. Put multi-word terms in quotes.

=item -i --ignore

Optional. Specify a factor type to ignore when creating assay groups.

=item -o --outdir

Optional. Specify a directory to write the XML configuration file. Default is
current working directory.

=item -a --assay-groups-only

Optional. Differential experiments only. Passing this option means the script
will only write assay groups to XML file, and will not try to create any
contrasts.

=item -d --debug

Optional. Log debugging messages.

=item -h --help

Print a helpful message. 

=back

=head1 AUTHOR

Expression Atlas team <arrayexpress-atlas@ebi.ac.uk>

=cut



use strict;
use warnings;
use 5.10.0;

use Pod::Usage;
use Getopt::Long;
use Cwd qw();
use DateTime;
use Log::Log4perl;
use Log::Log4perl::Level;
use File::Spec;
use Data::Dumper;

use Atlas::Common qw(
    create_non_strict_magetab4atlas
	create_atlas_site_config
);

use Atlas::AtlasConfig::Setup qw(
	create_factor_configs
	check_magetab4atlas_for_config
	create_atlas_experiment_type
);
use Atlas::AtlasConfig::ExperimentConfigFactory qw( create_experiment_config );


# Auto flush buffer.
$| = 1;

my $atlasSiteConfig = create_atlas_site_config;

# Parse command line arguments.
my $args = parse_args();

my $level = "INFO";

if( $args->{ "debug" } ) { $level = "DEBUG"; }

# Log4perl config.
my $logger_config = qq(
	log4perl.rootlogger			         = $level, LOG1, SCREEN
	log4perl.appender.SCREEN             = Log::Log4perl::Appender::Screen
	log4perl.appender.SCREEN.stderr      = 0
	log4perl.appender.SCREEN.layout      = Log::Log4perl::Layout::PatternLayout
	log4perl.appender.SCREEN.layout.ConversionPattern = %-5p - %m%n
	log4perl.appender.LOG1             = Log::Log4perl::Appender::File
	log4perl.appender.LOG1.filename    = sub { get_log_file_name }
	log4perl.appender.LOG1.header_text = sub { get_log_file_header }
	log4perl.appender.LOG1.mode        = append
	log4perl.appender.LOG1.layout      = Log::Log4perl::Layout::PatternLayout
	log4perl.appender.LOG1.layout.ConversionPattern = %-5p - %m%n
);

# Initialise logger.
Log::Log4perl::init(\$logger_config);
my $logger = Log::Log4perl::get_logger;

$logger->debug( "Debugging mode ON." );

my $atlasProdDir = $ENV{ "ATLAS_PROD" };

# Get name of file with reference factor values and factor types to ignore.
my $referencesIgnoreFile = File::Spec->catfile( $atlasProdDir, $atlasSiteConfig->get_references_ignore_file );

# Create hashes for reference factor values to use in contrasts, and factor
# types to ignore when creating assay groups.
$logger->info("Reading config for reference factor values and factor types to ignore from $referencesIgnoreFile");
my ($referenceFactorValues, $ignoreFactorTypes) = create_factor_configs($referencesIgnoreFile);

# If a manual reference factor value was specified, add it to the references hash.
if($args->{ "reference_value" }) {
	$logger->info("Using temporary reference value \"", $args->{ "reference_value" }, "\"");
	$referenceFactorValues->{ $args->{ "reference_value" } } = 1;
}

# If a manual factor type to ignore was specified, add it to the ignore hash.
if($args->{ "ignore_factor" }) {
	$logger->info("Temporarily ignoring factor type \"", $args->{ "ignore_factor" }, "\"");
	$ignoreFactorTypes->{ $args->{ "ignore_factor" } } = 1;
}

$logger->debug( "Parsing MAGE-TAB..." );

# Read in the MAGE-TAB.
my $magetab4atlas = create_non_strict_magetab4atlas( $args->{ "experiment_accession" } );

# Make sure the Magetab4Atlas object contains the appropriate assays.
$magetab4atlas = check_magetab4atlas_for_config($args, $ignoreFactorTypes, $magetab4atlas);

$logger->debug( "Successfully parsed MAGE-TAB" );

$logger->debug( Dumper( $magetab4atlas ) );

unless( $magetab4atlas->has_assays ) {
	$logger->logdie( "No assays were detected during MAGE-TAB parsing, cannot continue." );
}

# Create the XML config experiment type.
my $atlasExperimentType = create_atlas_experiment_type($magetab4atlas, $args->{ "analysis_type" });
$logger->info("Experiment type is $atlasExperimentType\n");

# Temporary: don't allow 2-colour experiments.
if( $atlasExperimentType =~ /2colour/ ) {
    $logger->logdie( "2-colour array experiments are not currently allowed." );
}

# Create the experiment config.
my $experimentConfig = create_experiment_config(
	$magetab4atlas, 
	$atlasExperimentType, 
	$args->{ "experiment_accession" }, 
	$referenceFactorValues,
	$args->{ "assay_groups_only" }
);

$experimentConfig->write_xml( $args->{ "output_directory" }, $args->{ "assay_groups_only" } );



###############
# Subroutines #
###############

# parse_args
# 	- Read arguments from command line.
sub parse_args {
	# Add the arguments to this hash with relevant options.
	my %args;
	# Set this if -h passed.
	my $want_help;
	
	# Possible analysis types.
	my @allowed_analysis_types = qw(
		baseline
		differential
	);
	# Possible library layouts.
	my @allowed_library_layouts = qw(
		single
		paired
	);	
	
	GetOptions(
		"h|help"			=> \$want_help,
		"e|experiment=s"	=> \$args{ "experiment_accession" },
		"t|type=s"			=> \$args{ "analysis_type" },	# baseline or differential
		"l|library=s"		=> \$args{ "library_layout" },	# paired or single
		"r|reference=s"		=> \$args{ "reference_value" },	# new reference factor value
		"i|ignore=s"		=> \$args{ "ignore_factor" },	# factor type to ignore
		"o|outdir=s"		=> \$args{ "output_directory" },	# dir for XML file
		"a|assay-groups-only"	=> \$args{ "assay_groups_only" },	# only create assay groups for differential (no contrasts)
		"d|debug"			=> \$args{ "debug" },
	);

	if($want_help) {
		pod2usage(
			-exitval => 255,
			-output => \*STDOUT,
			-verbose => 1
		);
	}

	# We must have experiment accession and type in order to do anything.
	unless($args{ "experiment_accession" } && $args{ "analysis_type" }) { 
		pod2usage(
			-message => "You must specify an experiment accession and an analysis type.\n",
			-exitval => 255,
			-output => \*STDOUT,
			-verbose => 1,
		);
	}

	# Check that accession is in the right format.
	unless($args{ "experiment_accession" } =~ /^E-\w{4}-\d+$/) {
		pod2usage(
			-message => "\"". $args{ "experiment_accession" }. "\" does not look like an ArrayExpress experiment accession.\n",
			-exitval => 255,
			-output => \*STDOUT,
			-verbose => 1,
		);

	}

	# Check that type is one of "baseline" or "differential".
	unless(grep $_ eq $args{ "analysis_type" }, @allowed_analysis_types) {
		pod2usage(
			-message => "\"". $args{ "analysis_type" }. "\" is not an allowed analysis type.\n",
			-exitval => 255,
			-output => \*STDOUT,
			-verbose => 1,
		);
	}

	# If we've been passed a library layout, check it's either "paired" or "single".
	if($args{ "library_layout" }) {
		unless(grep $_ eq $args{ "library_layout" }, @allowed_library_layouts) {
			pod2usage(
				-message => "\"". $args{ "library_layout" }. "\" is not an allowed library layout.\n",
				-exitval => 255,
				-output => \*STDOUT,
				-verbose => 1,
			);
		}
	}
	
	# If both "-t baseline" and "-r <value>" were passed, this doesn't make
	# sense -- don't need a reference for a baseline experiment. Die.
	if($args{ "analysis_type" } eq "baseline" && $args{ "reference_value" }) {
		pod2usage(
			-message => "Cannot use reference factor values in baseline experiments.\n",
			-exitval => 255,
			-output => \*STDOUT,
			-verbose => 1,
		);
	}

	# If "-a" is passed as well as "-t baseline", this is pointless, warn that
	# it doesn't make a difference but continue.
	if( $args{ "analysis_type" } eq "baseline" && $args{ "assay_groups_only" } ) {
		print "WARN  - assay groups only option makes no sense for a baseline experiment, ignoring.\n";
	}

	# If no output directory was specified, log that we will print to current working directory.
	unless($args{ "output_directory" }) {
		print "WARN  - No output directory specifed, will write XML configuration in ", Cwd::cwd(), "\n";
		$args{ "output_directory" } = Cwd::cwd();
	}

	# If one was specified, check that it's writable and die if not.
	unless(-w $args{ "output_directory" }) {
		pod2usage(
			-message => $args{ "output_directory" }. " is not writable or does not exist.\n",
			-exitval => 255,
			-output => \*STDOUT,
			-verbose => 1,
		);
	}

	return \%args;
}


sub get_log_file_name {
	
	my $logFileBaseName = "atlas_configuration_generation_" . $args->{ "experiment_accession" } . ".log";

	my $logFileName = File::Spec->catdir( $args->{ "output_directory" }, $logFileBaseName );

	# Delete the old one if it's there.
	if(-e $logFileName) {
		`rm $logFileName`;
	}


	return $logFileName;
}


sub get_log_file_header {

	my $headerText = "Atlas config generation log for " 
		. $args->{ "experiment_accession" } 
		. " created at " 
		. DateTime->now;
	
	$headerText .= "\n" . ("-" x 80) . "\n\n";
	
	return $headerText;
}

