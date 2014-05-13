#!/usr/bin/perl
#

use strict;
use warnings;

#######################################################################
# This here for testing.
use lib '/ebi/microarray/home/mkeays/Atlas/git/atlasprod/perl_modules';
#######################################################################

use AtlasConfig::Setup qw(
	create_factor_configs
	create_magetab4atlas
	create_atlas_experiment_type
);
use AtlasConfig::ExperimentConfigFactory qw( create_experiment_config );
use Getopt::Long;
use Cwd qw();
use Log::Log4perl;
use Log::Log4perl::Level;

# Auto flush buffer.
$| = 1;

# Log4perl config.
# We need to have two separate loggers because the info coming from commandline
# includes accession, which we want to use in the log file name, but we can't
# write to the file if e.g. no accession is provided, so we need to send that
# to the screen instead.
# Config for screen logger to print to screen...
my $screen_logger_config = q(
	log4perl.category.ATLASCONFIG_SCREEN = INFO, SCREEN
	log4perl.appender.SCREEN             = Log::Log4perl::Appender::Screen
	log4perl.appender.SCREEN.stderr      = 0
	log4perl.appender.SCREEN.layout      = Log::Log4perl::Layout::PatternLayout
	log4perl.appender.SCREEN.layout.ConversionPattern = %-5p - %m%n
);

# Config for file logger to print to a file.
my $file_logger_config = q(
	log4perl.category.ATLASCONFIG_FILE = INFO, LOG1
	log4perl.appender.LOG1             = Log::Log4perl::Appender::File
	log4perl.appender.LOG1.filename    = sub { get_log_file_name }
	log4perl.appender.LOG1.mode        = append
	log4perl.appender.LOG1.layout      = Log::Log4perl::Layout::PatternLayout
	log4perl.appender.LOG1.layout.ConversionPattern = %-5p - %m%n
);

# Initialise screen logger.
Log::Log4perl::init(\$screen_logger_config);
my $screen_logger = Log::Log4perl::get_logger("ATLASCONFIG_SCREEN");


#--------------------------------------------------
# use Log::Log4perl qw(:easy);
#-------------------------------------------------- 

#--------------------------------------------------
# Log::Log4perl->easy_init( { level => $INFO, layout => '%-5p - %m%n' } );
#-------------------------------------------------- 

# Parse command line arguments.
my $args = &parse_args();

# Now we have commandline args containing accession, initialise file logger as
# well.
#--------------------------------------------------
# Log::Log4perl::init(\$file_logger_config);
# my $file_logger = Log::Log4perl::get_logger("ATLASCONFIG_FILE");
#-------------------------------------------------- 

# Turn on debugging if required.
if($args->{ "debug" }) {
	#--------------------------------------------------
	# Log::Log4perl->easy_init( { level => $DEBUG, layout => '%-5p - %m%n' } );
	#-------------------------------------------------- 
	
	$screen_logger->level($DEBUG);

	$screen_logger->debug("Debugging mode ON");
}

# Hardcoding path to references/ignore file but FIXME.
# FIXME: Change back to $ATLAS_PROD one for production.
my $referencesIgnoreFile = "/ebi/microarray/home/mkeays/Atlas/jira/GRAMENE/gramene-62/reference_assay_group_factor_values.xml";

# Create hashes for reference factor values to use in contrasts, and factor
# types to ignore when creating assay groups.
$screen_logger->info("Reading config for reference factor values and factor types to ignore from $referencesIgnoreFile");
my ($referenceFactorValues, $ignoreFactorTypes) = create_factor_configs($referencesIgnoreFile);

if($args->{ "reference_value" }) {
	
	$screen_logger->info("Using temporary reference value \"", $args->{ "reference_value" }, "\"");

	$referenceFactorValues->{ $args->{ "reference_value" } } = 1;
}

# Get a Magetab4Atlas object containing the appropriate assays.
my $magetab4atlas = create_magetab4atlas($args, $ignoreFactorTypes);

# Create the XML config experiment type.
my $atlasExperimentType = create_atlas_experiment_type($magetab4atlas, $args->{ "analysis_type" });
$screen_logger->info("Experiment type is $atlasExperimentType\n");

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

	-x \"noreplicates\"
		Force configuration generation for experiments with no replicates.
		Don't do this unless you have to!
	
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
		"x|noreplicates"	=> \$args{ "no_replicates" },	# force allow no replicates
		"d|debug"			=> \$args{ "debug" },
	);

	if($want_help) {
		print $usage;
		exit;
	}

	# We must have experiment accession and type in order to do anything.
	unless($args{ "experiment_accession" } && $args{ "analysis_type" }) { 
		$screen_logger->error("Please specify \"-e <experiment accession> -t <baseline | differential>\""); 
	}

	# Check that accession is in the right format.
	unless($args{ "experiment_accession" } =~ /^E-\w{4}-\d+$/) {
		$screen_logger->logdie("\"", $args{ "experiment_accession" }, "\" does not look like an ArrayExpress experiment accession.");
	}

	# Check that type is one of "baseline" or "differential".
	unless(grep $_ eq $args{ "analysis_type" }, @allowed_analysis_types) {
		$screen_logger->logdie("\"", $args{ "analysis_type" }, "\" is not an allowed experiment type.");
	}

	# If we've been passed a library layout, check it's either "paired" or "single".
	if($args{ "library_layout" }) {
		unless(grep $_ eq $args{ "library_layout" }, @allowed_library_layouts) {
			$screen_logger->logdie("\"", $args{ "library_layout" }, "\" is not an allowed library layout.");
		}
	}
	
	# If both "-t baseline" and "-r <value>" were passed, this doesn't make
	# sense -- don't need a reference for a baseline experiment. Die.
	if($args{ "analysis_type" } eq "baseline" && $args{ "reference_value" }) {
		$screen_logger->logdie("Cannot use reference factor values in baseline experiments.");
	}

	# If no output directory was specified, log that we will print to current working directory.
	unless($args{ "output_directory" }) {
		$screen_logger->warn("No output directory specifed, will write XML configuration in ", Cwd::cwd() );
		$args{ "output_directory" } = Cwd::cwd();
	}
	# If one was specified, check that it's writable and die if not.
	unless(-w $args{ "output_directory" }) {
		$screen_logger->logdie($args{ "output_directory" }, " is not writable or does not exist.");
	}

	# If "noreplicates" is turned on, warn.
	if($args{ "no_replicates" }) {
		$screen_logger->warn("Allowing factor values with no biological replicates.");
	}
	
	return \%args;
}


sub get_log_file_name {
	my ($exptAccession) = @_;
	
	my $logFileName = "atlas_configuration_generation_".$exptAccession.".log";

	return $logFileName;
}











