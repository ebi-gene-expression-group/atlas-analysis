#!/usr/bin/perl
#

use strict;
use warnings;

#######################################################################
# This here for testing.
use lib '/ebi/microarray/home/mkeays/Atlas/git/atlasprod/perl_modules';
#######################################################################

use Getopt::Long;
use Cwd qw();
use DateTime;
use Log::Log4perl;
use Log::Log4perl::Level;

use AtlasConfig::Setup qw(
	create_factor_configs
	create_magetab4atlas
	create_atlas_experiment_type
);
use AtlasConfig::ExperimentConfigFactory qw( create_experiment_config );


# Auto flush buffer.
$| = 1;


# Parse command line arguments.
my $args = &parse_args();

# Log4perl config.
my $logger_config = q(
	log4perl.category.ATLASCONFIG_LOGGER          = INFO, LOG1, SCREEN
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
my $logger = Log::Log4perl::get_logger("ATLASCONFIG_LOGGER");

# Turn on debugging if required.
if($args->{ "debug" }) {
	$logger->level($DEBUG);
	$logger->debug("Debugging mode ON");
}

# Hardcoding path to references/ignore file but FIXME.
# FIXME: Change back to $ATLAS_PROD one for production.
my $referencesIgnoreFile = "/ebi/microarray/home/mkeays/Atlas/jira/GRAMENE/gramene-62/reference_assay_group_factor_values.xml";

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

# Get a Magetab4Atlas object containing the appropriate assays.
my $magetab4atlas = create_magetab4atlas($args, $ignoreFactorTypes);

# Create the XML config experiment type.
my $atlasExperimentType = create_atlas_experiment_type($magetab4atlas, $args->{ "analysis_type" });
$logger->info("Experiment type is $atlasExperimentType\n");

# Create the experiment config.
my $experimentConfig = create_experiment_config(
	$magetab4atlas, 
	$atlasExperimentType, 
	$args->{ "experiment_accession" }, 
	$referenceFactorValues
);

$experimentConfig->write_xml( $args->{ "output_directory" } );



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
	
	# A help message.
	my $usage = "
Usage:\n
	gxa_generateConfigurationForExperiment.pl -e <experiment accession> -t <baseline | differential> 

Options:
	-h \"help\"
		Print this help message.

	-e \"experiment\"
		ArrayExpress experiment accession e.g. \'E-MTAB-1066\'.

	-t \"type\"
		Type of analysis -- one of \'baseline\' or \'differential\'.

	-l \"library\"
		Specify the type of RNA-seq libraries to retrieve. One of
		\'paired\' or \'single\'.

	-r \"reference\"
		Differential experiments only. Specify a value to use as the
		reference factor value. E.g. \'-r \"untreated cells\"\'.

	-i \"ignore\"
		Specify a factor type to ignore. E.g. \'-i individual\'.

	-o \"outdir\"
		Specify the directory to write the XML configuration to. Default is
		current working directory.
	
	-d \"debug\"
		Print debugging messages.

";

	GetOptions(
		"h|help"			=> \$want_help,
		"e|experiment=s"	=> \$args{ "experiment_accession" },
		"t|type=s"			=> \$args{ "analysis_type" },	# baseline or differential
		"l|library=s"		=> \$args{ "library_layout" },	# paired or single
		"r|reference=s"		=> \$args{ "reference_value" },	# new reference factor value
		"i|ignore=s"		=> \$args{ "ignore_factor" },	# factor type to ignore
		"o|outdir=s"		=> \$args{ "output_directory" },	# dir for XML file
		"d|debug"			=> \$args{ "debug" },
	);

	if($want_help) {
		print $usage;
		exit;
	}

	# We must have experiment accession and type in order to do anything.
	unless($args{ "experiment_accession" } && $args{ "analysis_type" }) { 
		die "Please specify \"-e <experiment accession> -t <baseline | differential>\"\n"; 
	}

	# Check that accession is in the right format.
	unless($args{ "experiment_accession" } =~ /^E-\w{4}-\d+$/) {
		die "\"", $args{ "experiment_accession" }, "\" does not look like an ArrayExpress experiment accession.\n";
	}

	# Check that type is one of "baseline" or "differential".
	unless(grep $_ eq $args{ "analysis_type" }, @allowed_analysis_types) {
		die "\"", $args{ "analysis_type" }, "\" is not an allowed experiment type.\n";
	}

	# If we've been passed a library layout, check it's either "paired" or "single".
	if($args{ "library_layout" }) {
		unless(grep $_ eq $args{ "library_layout" }, @allowed_library_layouts) {
			die "\"", $args{ "library_layout" }, "\" is not an allowed library layout.\n";
		}
	}
	
	# If both "-t baseline" and "-r <value>" were passed, this doesn't make
	# sense -- don't need a reference for a baseline experiment. Die.
	if($args{ "analysis_type" } eq "baseline" && $args{ "reference_value" }) {
		die "Cannot use reference factor values in baseline experiments.\n";
	}

	# If no output directory was specified, log that we will print to current working directory.
	unless($args{ "output_directory" }) {
		print "WARN  - No output directory specifed, will write XML configuration in ", Cwd::cwd(), "\n";
		$args{ "output_directory" } = Cwd::cwd();
	}
	# If one was specified, check that it's writable and die if not.
	unless(-w $args{ "output_directory" }) {
		die $args{ "output_directory" }, " is not writable or does not exist.\n";
	}

	return \%args;
}


sub get_log_file_name {
	
	my $logFileName = $args->{ "output_directory" } . "atlas_configuration_generation_" . $args->{ "experiment_accession" } . ".log";

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








