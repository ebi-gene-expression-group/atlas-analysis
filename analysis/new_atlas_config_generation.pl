#!/usr/bin/perl
#

use strict;
use warnings;

#######################################################################
# This here for testing.
use lib '/ebi/microarray/home/mkeays/Atlas/git/atlasprod/perl_modules';
#######################################################################

use AtlasConfig::Setup qw(
	create_magetab4atlas
	create_atlas_experiment_type
);
use AtlasConfig::ExperimentConfigFactory qw( create_experiment_config );
use Getopt::Long;
use Cwd qw();
use Log::Log4perl qw(:easy);

Log::Log4perl->easy_init( { level => $INFO, layout => '%-5p - %m%n' } );

# Auto flush buffer.
$| = 1;

# Parse command line arguments.
my $args = &parse_args();

# Hardcoding path to references/ignore file but FIXME.
# FIXME: Change back to $ATLAS_PROD one for production.
my $referencesIgnoreFile = "/ebi/microarray/home/mkeays/Atlas/jira/GRAMENE/gramene-62/reference_assay_group_factor_values.txt";



# Get a Magetab4Atlas object containing the appropriate assays.
my $magetab4atlas = create_magetab4atlas($args);

# Create the XML config experiment type.
my $atlasExperimentType = create_atlas_experiment_type($magetab4atlas, $args->{ "analysis_type" });
INFO "Experiment type is $atlasExperimentType\n";

# Create the AtlasConfig::ExperimentConfig
my $experimentConfig = create_experiment_config($magetab4atlas, $atlasExperimentType, $args->{ "experiment_accession" });

$experimentConfig->write_xml($args->{ "output_directory" });


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
	);

	if($want_help) {
		print $usage;
		exit;
	}

	# We must have experiment accession and type in order to do anything.
	unless($args{ "experiment_accession" } && $args{ "analysis_type" }) { 
		ERROR "Please specify \"-e <experiment accession> -t <baseline | differential>\""; 
	}

	# Check that accession is in the right format.
	unless($args{ "experiment_accession" } =~ /^E-\w{4}-\d+$/) {
		LOGDIE "\"", $args{ "experiment_accession" }, "\" does not look like an ArrayExpress experiment accession.";
	}

	# Check that type is one of "baseline" or "differential".
	unless(grep $_ eq $args{ "analysis_type" }, @allowed_analysis_types) {
		LOGDIE "\"", $args{ "analysis_type" }, "\" is not an allowed experiment type.";
	}

	# If we've been passed a library layout, check it's either "paired" or "single".
	if($args{ "library_layout" }) {
		unless(grep $_ eq $args{ "library_layout" }, @allowed_library_layouts) {
			LOGDIE "\"", $args{ "library_layout" }, "\" is not an allowed library layout.";
		}
	}
	
	# If both "-t baseline" and "-r <value>" were passed, this doesn't make
	# sense -- don't need a reference for a baseline experiment. Die.
	if($args{ "analysis_type" } eq "baseline" && $args{ "reference_value" }) {
		LOGDIE "Cannot use reference factor values in baseline experiments.";
	}

	# If no output directory was specified, log that we will print to current working directory.
	unless($args{ "output_directory" }) {
		WARN "No output directory specifed, will write XML configuration in ", Cwd::cwd();
		$args{ "output_directory" } = Cwd::cwd();
	}

	# If "noreplicates" is turned on, warn.
	if($args{ "no_replicates" }) {
		WARN "Allowing factor values with no biological replicates.";
	}
	
	return \%args;
}

