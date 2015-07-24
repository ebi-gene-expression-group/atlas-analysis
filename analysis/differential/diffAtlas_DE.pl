#!/usr/bin/env perl
#
# diffAtlas_DE.pl
# 	- Get results of differential expression for either microarray or RNA-seq data.
# 	- For microarray data, this entails finding contrast definitions in contrasts XML file, running
# 	limma for each contrast, creating an MvA plot for each contrast, and then
# 	summarizing statistics results for all contrasts into a single
# 	tab-delimited text file.
# 	- For RNA-seq, differential expression stats for each contrast have already
# 	been calculated using DESeq via iRAP. For this data, find contrast
# 	definitions in iRAP config and display names in Atlas contrasts XML file,
# 	read DESeq results for each contrast, create an MvA plot for each contrast,
# 	and then summarize statistics results for all contrasts into a single
# 	tab-delimited text file.

=pod

=head1 NAME

diffAtlas_DE.pl - do differential expression analysis for Expression Atlas.

=head1 SYNOPSIS

diffAtlas_DE.pl --experiment E-MTAB-1066 --directory path/to/processing/directory

diffAtlas_DE.pl --experiment E-GEOD-38400 --irapconf E-MTAB-38400.conf

=head1 DESCRIPTION

This script currently either (a) runs differential expression analysis via
limma for microarray experiments, or (b) reads in iRAP's DESeq results files
and compiles them for an RNA-seq experiment.

=head1 OPTIONS

=over 2

=item -e --experiment

Required. ArrayExpress accession of experiment.

=item -d --directory

Required. Path to Atlas processing directory containin XML config and
normalized data files (for microarray).

=item -i --irapconf

Required for RNA-seq experiments. Full path to iRAP config file.

=item -h --help

Optional. Print a helpful message.

=back

=head1 AUTHOR

Expression Atlas team <arrayexpress-atlas@ebi.ac.uk

=cut


use strict;
use warnings;
use 5.10.0;

use Cwd;
use File::Spec;
use File::Basename;
use Getopt::Long;
use IPC::Cmd qw( can_run );
use Log::Log4perl;
use Pod::Usage;
use Scalar::Util qw(looks_like_number);

use Atlas::AtlasConfig::Reader qw( parseAtlasConfig );
use Atlas::AtlasConfig::Common qw( map_technical_replicate_ids_to_assays );


# Flush after every print not every newline
$| = 1;

my $logger_config = q(
	log4perl.rootlogger			         = INFO, SCREEN
	log4perl.appender.SCREEN             = Log::Log4perl::Appender::Screen
	log4perl.appender.SCREEN.stderr      = 0
	log4perl.appender.SCREEN.layout      = Log::Log4perl::Layout::PatternLayout
	log4perl.appender.SCREEN.layout.ConversionPattern = %-5p - %m%n
);

# Initialise logger.
Log::Log4perl::init(\$logger_config);
my $logger = Log::Log4perl::get_logger;


# assumes R scripts directory is in PATH
my $mvaScript = "diffAtlas_mvaPlot.R";
my $limmaScript = "diffAtlas_DE_limma.R";
my $deseqScript = "diffAtlas_DE_deseq.R";

# Check that R is installed.
unless( can_run( "Rscript" ) ) {
    $logger->logdie( "R not found. Please ensure it is installed and you can run it." );
}

# Check the R scripts can be run.
foreach my $Rscript ( $mvaScript, $limmaScript, $deseqScript ) {

    unless( can_run( $Rscript ) ) {
        $logger->logdie( "$Rscript not found. Please ensure it is installed and you can run it." );
    }
}

# Get commandline arguments
my $args = parse_args();

my $atlasXML = File::Spec->catfile( 
	$args->{ "processing_directory" }, 
	$args->{ "experiment_accession" } . "-configuration.xml" 
);

# Read XML config into Atlas::AtlasConfig::ExperimentConfig object.
my $experimentConfig = parseAtlasConfig( $atlasXML );

my $atlasExperimentType = $experimentConfig->get_atlas_experiment_type;

# If this is an RNA-seq experiment, make sure we got the path to the iRAP config file.
if( $atlasExperimentType =~ "rnaseq" ) {
	unless( $args->{ "irap_config" } ) {
		$logger->logdie( 
			$args->{ "experiment_accession" }, 
			" is an RNA-seq experiment but no iRAP config file was specified. Cannot continue."
		);
	}
}

# Make sure we haven't been given a baseline experiment.
if( $atlasExperimentType =~ /baseline/ ) {
	$logger->logdie( "Experiments with type \"", $experimentConfig->get_atlas_experiment_type, "\" cannot be processed as differential." );
}


########################################
# If everything is OK, begin processing!
# For microarray experiments, run limma analysis, write results, and create MvA
# plots.
if( $atlasExperimentType =~ /array/ ) {
	
	run_microarray_differential_expression( 
		$limmaScript,
		$args->{ "experiment_accession" }, 
		$experimentConfig, 
		$args->{ "processing_directory" }, 
	);
}
# For RNA-seq experiments, take results from iRAP's DESeq results files, write
# them into one file, and create MvA plots.
elsif( $atlasExperimentType =~ /rnaseq/ ) {
	
	# Read the iRAP config file to get the path to the raw counts matrix.
	my $countsMatrixFilePath = get_counts_matrix_path( $args->{ "irap_config" }, $args->{ "experiment_accession" } );
	
	run_rnaseq_differential_expression(
		$deseqScript,
		$args->{ "experiment_accession" },
		$experimentConfig,
		$args->{ "processing_directory" },
		$countsMatrixFilePath
	);

	# Get the results, write them out, make MvA plots.
	#get_rnaseq_results( $contrastIDsToDESeqFiles, $experimentConfig, $args->{ "experiment_accession" }, $args->{ "processing_directory" } );
}
# If we're passed an experiment type we don't recognise, die.
else {
	$logger->logdie( "Don't know what to do with experiment of type \"", $atlasExperimentType, "\"\n" );
}


###############
# subroutines #
###############

# init
# 	Get commandline arguments using Getopt::Long.
# 		--atlasxml <Atlas contrasts XML file> (required)
# 		--irapconf <iRAP config file> (for RNA-seq)
sub parse_args {
	
	my %args;

	my $want_help;

	GetOptions(
		"h|help"			=> \$want_help,
		"e|experiment=s"	=> \$args{ "experiment_accession" }, 
		"d|directory=s"		=> \$args{ "processing_directory" },
		"i|irapconf=s"		=> \$args{ "irap_config" },
	);

	# Print help if requested.
	if( $want_help ) {
		pod2usage(
			-exitval	=> 255,
			-output		=> \*STDOUT,
			-verbose	=> 1
		);
	}

	# Check that we have the minimum requirements.
	unless( $args{ "experiment_accession" } && $args{ "processing_directory" } ) {
		pod2usage(
			-message	=> "You must specify an experiment accession and a processing directory.\n",
			-exitval	=> 255,
			-output		=> \*STDOUT,
			-verbose	=> 1
		);
	}

	# Check that the experiment accession is a valid one.
	unless( $args{ "experiment_accession" } =~ /^E-\w{4}-\d+$/ ) {
		pod2usage(
			-message	=> $args{ "experiment_accession" } . " does not look like an ArrayExpress experiment accession.\n",
			-exitval	=> 255,
			-output		=> \*STDOUT.
			-verbose	=> 1
		);
	}

	return \%args;
}


# run_microarray_differential_expression
# 	- Run differential expression analysis using limma script.
# 	- For each array design, combine results for all contrasts and write
# 	adjusted p-values, moderated t-statistics, and log2(fold-change)s to one
# 	file.
# 	- Create MvA plots for each contrast using MvA plotting script.
sub run_microarray_differential_expression {

	my ( $limmaScript, $expAcc, $experimentConfig, $atlasProcessingDirectory ) = @_;

	$logger->info( "Running differential expression analysis in R..." );
	
    my $tempDir = File::Spec->catdir( $ENV{ "HOME" }, "tmp" );

	unless( -d $tempDir ) {
		$logger->logdie( "No $tempDir found. Please create it and try again." );
	}
	
	# Run R script.
	my $Routput = `$limmaScript $expAcc $atlasProcessingDirectory 2>&1`;

	# Check R output for errors.
	if( $? ) {

		# Can't continue without results from limma so may as well quit.
		$logger->logdie( "Problems during differential expression analysis:\n\n$Routput" );

	} else {
		$logger->info( "Differential expression analysis successful" );
	}
			
	# Get the differential expression results into the hash of all analytics results.
	# Files are read out of the users ~/tmp directory where the R script wrote
	# them to.
	my $analyticsDEResults = read_limma_results( $expAcc );

	# Map contrast IDs to contrast names.
	my $contrastIDs2names = map_contrast_ids_to_names( $experimentConfig );
	my $contrastIDs2arrayDesigns = map_contrast_ids_to_arraydesigns( $experimentConfig );

	# Get the names of the MvA plot data files.
	my @plotDataFiles = glob( "$tempDir/$expAcc.g*_g*.plotdata.tsv" );

    unless( @plotDataFiles ) {
        $logger->logdie( "No MvA plot data files found in $tempDir , please check." );
    }

	foreach my $plotDataFile ( @plotDataFiles ) {

        $logger->info( "Reading file $plotDataFile ..." );

		( my $contrastID = basename( $plotDataFile ) ) =~ s/.*\.(g\d+_g\d+)\.plotdata.*/$1/;

		my $contrastName = $contrastIDs2names->{ $contrastID };
		my $arrayDesignAccession = $contrastIDs2arrayDesigns->{ $contrastID };

		# Filename for MvA plot.
		my $plotFile = File::Spec->catfile( $atlasProcessingDirectory, $expAcc."_".$arrayDesignAccession."-".$contrastID."-mvaPlot.png" );

		# Create MvA.
		make_mva_plot( "microarray", $plotFile, $plotDataFile, $contrastName, $mvaScript );

        $logger->info( "Finished with file $plotDataFile" );
	}

	# Now we have results for all the contrasts in this analytics element. Write them to a file.
	write_results( $analyticsDEResults, $expAcc, $experimentConfig, $atlasProcessingDirectory );
	
}


sub map_contrast_ids_to_arraydesigns {

	my ( $experimentConfig ) = @_;

	my $allAnalytics = $experimentConfig->get_atlas_analytics;

	my $contrastIDs2arrayDesigns = {};

	foreach my $analytics ( @{ $allAnalytics } ) {

		my $contrasts = $analytics->get_atlas_contrasts;
		my $platform = $analytics->get_platform;

		foreach my $contrast ( @{ $contrasts } ) {

			my $contrastID = $contrast->get_contrast_id;

			$contrastIDs2arrayDesigns->{ $contrastID } = $platform;
		}
	}

	return $contrastIDs2arrayDesigns;
}


sub map_contrast_ids_to_names {

	my ( $experimentConfig ) = @_;

	my $allAnalytics = $experimentConfig->get_atlas_analytics;

	my $contrastIDs2names = {};

	foreach my $analytics ( @{ $allAnalytics } ) {

		my $contrasts = $analytics->get_atlas_contrasts;

		foreach my $contrast ( @{ $contrasts } ) {

			my $contrastID = $contrast->get_contrast_id;
			my $contrastName = $contrast->get_contrast_name;

			$contrastIDs2names->{ $contrastID } = $contrastName;
		}
	}

	return $contrastIDs2names;
}

# read_limma_results
# 	- parses output from limma script and adds them to a hash.
sub read_limma_results {

	my ( $expAcc ) = @_;

	my $analyticsDEResults = {};

	my $tempDir = File::Spec->catdir( $ENV{ "HOME" }, "tmp" );
	
	unless( -d $tempDir ) {
		$logger->logdie( "No $tempDir found. Please create it and try again." );
	}

	my @limmaResultsFiles = glob( "$tempDir/$expAcc.g*_g*.analytics.tsv" );

    unless( @limmaResultsFiles ) {
        $logger->logdie( "No limma results files found in $tempDir , please check." );
    }

	foreach my $limmaResultsFile ( @limmaResultsFiles ) {

        $logger->info( "Reading file $limmaResultsFile ..." );

		( my $contrastID = basename( $limmaResultsFile ) ) =~ s/.*\.(g\d+_g\d+)\.analytics.*/$1/;
	
		# Add results to hash.
		open( my $limmaResFH, "<$limmaResultsFile" ) or $logger->logdie( "Can't open file $limmaResultsFile: $!" );
		while(defined(my $line = <$limmaResFH>)) {
			
			# potentially don't need headers in the limma results file?
			unless($line =~ /^designElements/) {
				
				chomp $line;
				
				# Split on tabs
				my @lineSplit = split "\t", $line;
				
				# Design element (probeset) ID is the first element, take it off.
				my $designElement = shift @lineSplit;
				
				# Add array [p-value, t-statistic, logFoldChange] to hash for this
				# design element for this contrast.
				$analyticsDEResults->{ $designElement }->{ $contrastID } = \@lineSplit;
			}
		}
		close($limmaResFH);
        
        $logger->info( "Finished with $limmaResultsFile , deleting..." );

		# clean up
		`rm $limmaResultsFile`;

        $logger->info( "Done." );
	}

	return $analyticsDEResults;
}


# get_counts_matrix_path
#  - read iRAP config file to get filename for raw counts matrix.
sub get_counts_matrix_path {

	my ( $irapConfig, $expAcc ) = @_;

	# Directory of config file is the same as the beginning of path to
	# DESeq files.
	my $irapDir = ( fileparse( $irapConfig ) )[1];

	# Get experiment accession, mapper and quantification method
	# from config file to create full path to counts matrix file.
	my ( $irapExptAcc, $mapper, $qMethod );
	
	open( my $confFH, $irapConfig ) or die( "Can't open $irapConfig: $!\n" );

	$logger->info( "Reading iRAP config: $irapConfig..." );
	while( defined( my $line = <$confFH> ) ) {

		chomp $line;

		if( $line =~ /^name=(.*)/ ) { 
			$irapExptAcc = $1; 
		
			# Check experiment accession from iRAP config matches that in filename of
			# Atlas contrasts file.
			unless( $irapExptAcc eq $expAcc ) {
				$logger->logdie( "iRAP experiment accession is: $irapExptAcc -- does not match Atlas experiment accession ($expAcc)" );
			}
		}
		
		elsif( $line =~ /^mapper=(.*)/ ) { $mapper = $1; }
		elsif( $line =~ /^quant_method=(.*)/ ) { $qMethod = $1; }
	}
	close( $confFH );

	# Build path to raw counts matrix.
	my $countsMatrixFilePath = File::Spec->catfile( $irapDir, $expAcc, $mapper, $qMethod, "genes.raw." . $qMethod . ".tsv" );
	
	# Make sure the file exists and is readable.
	unless( -r $countsMatrixFilePath ) {
		$logger->logdie( 
			"Cannot read counts matrix for $expAcc from $countsMatrixFilePath\nPlease verify this file exists and is readable."
		);
	}
	
	# If all's OK, return the file path.
	return $countsMatrixFilePath;
}


# run_rnaseq_differential_expression
# 	- Run differential expression analysis using DESeq script.
# 	- Combine results for all contrasts and write adjusted p-values, moderated
# 	t-statistics, and log2(fold-change)s to one file.
# 	- Create MvA plots for each contrast using MvA plotting script.
sub run_rnaseq_differential_expression {

	my ( $deseqScript, $expAcc, $experimentConfig, $atlasProcessingDir, $countsMatrixFile ) = @_;

	$logger->info( "Running differential expression analysis in R..." );
	
    my $tempDir = File::Spec->catdir( $ENV{ "HOME" }, "tmp" );

	unless( -d $tempDir ) {
		$logger->logdie( "No $tempDir found. Please create it and try again." );
	}
	
	# Run R script.
	my $Routput = `$deseqScript $expAcc $atlasProcessingDir $countsMatrixFile 2>&1`;

	# Check R output for errors.
	if( $? ) {

		# Can't continue without results from limma so may as well quit.
		$logger->logdie( "Problems during differential expression analysis:\n\n$Routput" );

	} else {
		$logger->info( "Differential expression analysis successful" );
	}
			
	# Get the differential expression results into the hash of all analytics results.
	# Files are read out of the users ~/tmp directory where the R script wrote
	# them to.
	my $analyticsDEResults = read_deseq_results( $expAcc );

	# Map contrast IDs to contrast names.
	my $contrastIDs2names = map_contrast_ids_to_names( $experimentConfig );

	# Get the names of the MvA plot data files.
	my @plotDataFiles = glob( "$tempDir/$expAcc.g*_g*.plotdata.tsv" );

    unless( @plotDataFiles ) {
        $logger->logdie( "No MvA plot data files found in $tempDir , please check." );
    }

	foreach my $plotDataFile ( @plotDataFiles ) {

        $logger->info( "Reading file $plotDataFile ..." );

		( my $contrastID = basename( $plotDataFile ) ) =~ s/.*\.(g\d+_g\d+)\.plotdata.*/$1/;

		my $contrastName = $contrastIDs2names->{ $contrastID };

		# Filename for MvA plot.
		my $plotFile = File::Spec->catfile( $atlasProcessingDir, $expAcc."-".$contrastID."-mvaPlot.png" );

		# Create MvA.
		make_mva_plot( "rnaseq", $plotFile, $plotDataFile, $contrastName, $mvaScript );

        $logger->info( "Finished with file $plotDataFile" );
	}

	# Now we have results for all the contrasts in this analytics element. Write them to a file.
	write_results( $analyticsDEResults, $expAcc, $experimentConfig, $atlasProcessingDir );
	
}


sub read_deseq_results {

	my ( $expAcc ) = @_;
	
	my $analyticsDEResults = {};

	my $tempDir = File::Spec->catdir( $ENV{ "HOME" }, "tmp" );
	
	unless( -d $tempDir ) {
		$logger->logdie( "No $tempDir found. Please create it and try again." );
	}

	my @deseqResultsFiles = glob( "$tempDir/$expAcc.g*_g*.analytics.tsv" );

    unless( @deseqResultsFiles ) {
        $logger->logdie( "No DESeq results files found in $tempDir , please check." );
    }

	foreach my $deseqResultsFile ( @deseqResultsFiles ) {

        $logger->info( "Reading file $deseqResultsFile ..." );

		( my $contrastID = basename( $deseqResultsFile ) ) =~ s/.*\.(g\d+_g\d+)\.analytics.*/$1/;
	
		my ( $logfcIdx, $adjpvalIdx );
		
		# Add results to hash.
		open( my $deseqResFH, "<$deseqResultsFile" ) or $logger->logdie( "Can't open file $deseqResultsFile: $!" );
		while(defined(my $line = <$deseqResFH>)) {
			
			chomp $line;
			
			# Split on tabs
			my @lineSplit = split "\t", $line;

			# On first line, get column indices
			if($line =~ /^id\t/) {

				$logfcIdx = (grep { $lineSplit[$_] eq "log2FoldChange" } 0..$#lineSplit)[0];
				$adjpvalIdx = (grep { $lineSplit[$_] eq "padj" } 0..$#lineSplit)[0];
			}
			else {
				
				# Use indices found above to get values (gene ID is always the
				# very first element).
				my $geneID = $lineSplit[0];

				my $logFC = $lineSplit[$logfcIdx];
				my $adjPval = $lineSplit[$adjpvalIdx];
				# Ensure we get numbers for the things that are supposed to
				# be numbers (weird prob with DESeq output).
				# Allow "NA" adjusted p-values through.
				unless(looks_like_number($adjPval) || $adjPval eq "NA") {
					$logger->logdie( "\nDid not get numeric value for adjusted p-value:\nGene ID: $geneID\nadjusted p-value: $adjPval" );
				}
				# Allow "NA" log fold-changes through.
				unless(looks_like_number($logFC) || $logFC eq "NA") {
					$logger->logdie( "\nDid not get numeric value for log2FoldChange:\nGene ID: $geneID\nlog2FoldChange: $logFC" );
				}
				
				# Add to hash for file of all contrasts' results.
				$analyticsDEResults->{ $geneID }->{ $contrastID } = [ $adjPval, $logFC ];
			}
		}

		close($deseqResFH);
        
        $logger->info( "Finished with $deseqResultsFile , deleting..." );

		# clean up
		`rm $deseqResultsFile`;

        $logger->info( "Done." );
	}

	return $analyticsDEResults;
}


# get_rnaseq_results
# 	- For each contrast, parse the DESeq results file to get log2(fold-change)s
# 	and adjusted p-values.
# 	- Write all contrast results into one single file.
# 	- For each contrast, create MvA plots using MvA plotting script.
#--------------------------------------------------
# sub get_rnaseq_results {
# 
# 	my ( $contrastIDsToDESeqFiles, $experimentConfig, $expAcc, $atlasProcessingDirectory ) = @_;
# 
# 	# Get the RNA-seq analytics section from the experiment config.
# 	my $allAnalytics = $experimentConfig->get_atlas_analytics;
# 	
# 	# Empty variable to put RNA-seq analytics object in.
# 	my $rnaSeqAnalytics;
# 	
# 	# Go through the analytics elements, get the rnaseq one.
# 	foreach my $analytics ( @{ $allAnalytics } ) {
# 		
# 		# Get the contrasts from the rnaseq analytics element.
# 		if( $analytics->get_platform eq "rnaseq" ) { $rnaSeqAnalytics = $analytics; }
# 	}	
# 
# 	# Get the contrasts from the RNA-seq analytics object.
# 	my $rnaSeqContrasts = $rnaSeqAnalytics->get_atlas_contrasts;
# 
# 	# Make a hash for the contrasts mapping IDs to names.
# 	my $contrastIDsToNames = {};
# 	foreach my $contrast ( @{ $rnaSeqContrasts } ) {
# 
# 		my $contrastID = $contrast->get_contrast_id;
# 		my $contrastName = $contrast->get_contrast_name;
# 
# 		$contrastIDsToNames->{ $contrastID } = $contrastName;
# 	}
# 
# 	# Empty hash to store differential expression results before writing.
# 	my $rnaSeqDEResults = {};
# 
# 	$logger->info( "\nCollecting DESeq results..." );
# 	
# 	# Go through the contrasts...
# 	foreach my $contrastID (keys %{ $contrastIDsToDESeqFiles } ) {
# 
# 		$logger->info( "Collecting results for contrast: ", $contrastIDsToNames->{ $contrastID } );
# 		
# 		# Get the path to the DESeq results file.
# 		my $deseqFile = $contrastIDsToDESeqFiles->{ $contrastID };
# 
# 		# Need to get column indices of "baseMean", "log2FoldChange", and "padj" in DESeq results.
# 		my ($basemeanIdx, $logfcIdx, $adjpvalIdx);
# 
# 		open(DESEQRES, $deseqFile) or $logger->logdie( "Can't open $deseqFile: $!" );
# 
# 		# Temp file for MvA plot data.
# 		my $plotDataTempFile = "/tmp/plotData.$$.txt";
# 		
# 		open(PLOTDATA, ">$plotDataTempFile") or $logger->logdie( "Can't open file to write temporary plot data to: $!" );
# 		printf(PLOTDATA "ID\tavgExpr\tlogFC\tadjPval");
# 
# 		print "Reading DESeq results from $deseqFile...";
# 		while(defined(my $line = <DESEQRES>)) {
# 
# 			chomp $line;
# 			my @lineSplit = split "\t", $line;
# 			
# 			# On first line, get column indices
# 			if($line =~ /^id\t/) {
# 
# 				$basemeanIdx = (grep { $lineSplit[$_] eq "baseMean" } 0..$#lineSplit)[0];
# 				$logfcIdx = (grep { $lineSplit[$_] eq "log2FoldChange" } 0..$#lineSplit)[0];
# 				$adjpvalIdx = (grep { $lineSplit[$_] eq "padj" } 0..$#lineSplit)[0];
# 
# 				if($basemeanIdx == 0 || $logfcIdx == 0 || $adjpvalIdx == 0) {
# 					$logger->logdie( "Couldn't get column indices for all required columns." );
# 				}
# 			}
# 			else {
# 				
# 				# Use indices found above to get values (gene ID is always the
# 				# very first element).
# 				my $geneID = $lineSplit[0];
# 
# 				my $baseMean = $lineSplit[$basemeanIdx];
# 				my $logFC = $lineSplit[$logfcIdx];
# 				my $adjPval = $lineSplit[$adjpvalIdx];
# 				# Ensure we get numbers for the things that are supposed to
# 				# be numbers (weird prob with DESeq output).
# 				unless(looks_like_number($baseMean)) {
# 					$logger->logdie( "\nDid not get numeric value for baseMean:\nGene ID: $geneID\nbaseMean: $baseMean" );
# 				}
# 				# Allow "NA" adjusted p-values through.
# 				unless(looks_like_number($adjPval) || $adjPval eq "NA") {
# 					$logger->logdie( "\nDid not get numeric value for adjusted p-value:\nGene ID: $geneID\nadjusted p-value: $adjPval" );
# 				}
# 				# Allow "NA" log fold-changes through.
# 				unless(looks_like_number($logFC) || $logFC eq "NA") {
# 					$logger->logdie( "\nDid not get numeric value for log2FoldChange:\nGene ID: $geneID\nlog2FoldChange: $logFC" );
# 				}
# 				
# 				# Add to hash for file of all contrasts' results.
# 				$rnaSeqDEResults->{ $geneID }->{ $contrastID } = [ $adjPval, $logFC ];
# 				# Add to file for MvA plot.
# 				printf(PLOTDATA "\n$geneID\t$baseMean\t$logFC\t$adjPval");
# 			}
# 		}
# 		close(PLOTDATA);
# 		close(DESEQRES);
# 		$logger->info( "done" );
# 
# 		# Filename for MvA plot
# 		my $plotFile = "$expAcc-".$contrastID."-mvaPlot.png";
# 
# 		# Get the human-readable contrast name from the hash.
# 		my $atlasName = $contrastIDsToNames->{ $contrastID };
# 
# 		# Create MvA
# 		make_mva_plot("rnaseq", $plotFile, $plotDataTempFile, $atlasName, $mvaScript);
# 	}
# 	
# 	# Write the results to a file.
# 	write_results( $rnaSeqDEResults, $experimentConfig->get_experiment_accession, $experimentConfig, $atlasProcessingDirectory );
# }
#-------------------------------------------------- 


# make_mva_plot
# 	- Create MvA plot using MvA plot script.
sub make_mva_plot {

	my ( $platform, $plotFile, $plotDataTempFile, $contrastName, $mvaScript ) = @_;

	$logger->info( "Making MvA plot..." );
	# Create MvA plot with MvA plot script
	my $R_mvaOutput = `$mvaScript $plotDataTempFile \"$contrastName\" $plotFile $platform 2>&1`;
	
	# Check for errors (this doesn't catch warnings).
	if( $? ) {

		# Report the error but don't worry about dying as we can live
		# without the MvA plot.
		$logger->logdie( "Problems creating MvA plot:\n\n$R_mvaOutput" );

	} else {
		$logger->info( "done." );
	}
	# clean up
	`rm $plotDataTempFile`;
	`rm Rplots.pdf`;
}


# write_results
#	- write all results for a single platform (e.g. "A-AFFY-35" or "rnaseq") to
#	a tab-delimited text file.
sub write_results {

	my ( $analyticsDEResults, $expAcc, $experimentConfig, $atlasProcessingDirectory ) = @_;
	
	foreach my $analytics ( @{ $experimentConfig->get_atlas_analytics } ) {

		# Get the platform from the analytics element.
		my $platform = $analytics->get_platform;

		# Results file needs array design accession if it's microarray.
		my $resFile = ( $platform eq "rnaseq" ? $expAcc."-analytics.tsv.undecorated" : $expAcc."_".$platform."-analytics.tsv.undecorated" );

		$resFile = File::Spec->catfile( $atlasProcessingDirectory, $resFile );

		# Get all the contrast IDs for this analytics element.
		my $analyticsContrastIDs = [];
		
		foreach my $contrast ( @{ $analytics->get_atlas_contrasts } ) {

			push @{ $analyticsContrastIDs }, $contrast->get_contrast_id;
		}
		
		$logger->info( "Writing results to $resFile" );

		open( my $resFH, ">$resFile") or $logger->logdie( "Can't open $resFile: $!" );

		# Write column headers
		if($platform eq "rnaseq") { printf($resFH "Gene ID"); }
		else { printf($resFH "Design Element"); }
		foreach my $contrastID ( @{ $analyticsContrastIDs }) {
			printf($resFH "\t$contrastID.p-value");
			unless($platform eq "rnaseq") { printf($resFH "\t$contrastID.t-statistic"); }
			printf($resFH "\t$contrastID.log2foldchange");
		}

		# Write statistics
		foreach my $id (sort keys %{ $analyticsDEResults }) {

			printf($resFH "\n$id");

			# Use ordering of assay group pairs from %contrastHash so stats are in
			# the same order as the headers.
			foreach my $contrastID ( @{ $analyticsContrastIDs }) {

				my $statsString = "";
				if(exists($analyticsDEResults->{ $id }->{ $contrastID })) {
					$statsString = join "\t", @{ $analyticsDEResults->{ $id }->{ $contrastID }};
				} else {

					# Some genes are excluded if their counts are zero across all
					# assays. If there is more than one contrast in the experiment,
					# assays in one contrast may have non-zero counts while assays
					# in another contrast all have zero counts. In this case, we
					# have DESeq statistics results for one contrast and none for
					# the other. Where we don't have any results, we will put "NA"
					# values into the results file.
					if($platform eq "rnaseq") { $statsString = "NA\tNA"; }
					else { $statsString = "NA\tNA\tNA"; }
				}
				printf($resFH "\t$statsString");
			}
		}
		close( $resFH );
	}
}
			

