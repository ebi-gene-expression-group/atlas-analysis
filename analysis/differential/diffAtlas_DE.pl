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

use strict;
use warnings;

use File::Spec;
use File::Basename;
use Scalar::Util qw(looks_like_number);
use IPC::Cmd qw( can_run );
use Cwd;

use Atlas::AtlasConfig::Reader qw( parseAtlasConfig );
use Atlas::AtlasConfig::Common qw( map_technical_replicate_ids_to_assays );

# Flush after every print not every newline
$| = 1;

# assumes R scripts directory is in PATH
my $mvaScript = "diffAtlas_mvaPlot.R";
my $limmaScript = "diffAtlas_DE_limma.R";

# Check that R is installed.
unless( can_run( "R" ) ) {
    die( "R not found. Please ensure it is installed and you can run it.\n" );
}

# Check the R scripts can be run.
foreach my $Rscript ( $mvaScript, $limmaScript ) {

    unless( can_run( $Rscript ) ) {
        die( "$Rscript not found. Please ensure it is installed and you can run it.\n" );
    }
}

# Get commandline arguments
my ($atlasXML, $irapConfig) = init();


# Read XML config into Atlas::AtlasConfig::ExperimentConfig object.
my $experimentConfig = parseAtlasConfig( $atlasXML );

my $atlasExperimentType = $experimentConfig->get_atlas_experiment_type;

# Make sure we haven't been given a baseline experiment.
if( $atlasExperimentType =~ /baseline/ ) {
	die( "Experiments with type \"", $experimentConfig->get_atlas_experiment_type, "\" cannot be processed as differential.\n" );
}

# Get the experiment accession.
my $expAcc = $experimentConfig->get_experiment_accession;

print "Experiment accession: $expAcc\n";

# For microarray experiments, run limma analysis, write results, and create MvA
# plots.
if( $atlasExperimentType =~ /array/ ) {
	
	run_microarray_differential_expression( $expAcc, $experimentConfig );
}

# For RNA-seq experiments, take results from iRAP's DESeq results files, write
# them into one file, and create MvA plots.
elsif( $atlasExperimentType =~ /rnaseq/ ) {
	
	# Read the iRAP config file to get the paths to the DESeq results files for
	# each contrast.
	my $contrastIDsToDESeqFiles = read_irap_conf( $irapConfig, $expAcc );
	
	# Get the results, write them out, make MvA plots.
	get_rnaseq_results( $contrastIDsToDESeqFiles, $experimentConfig );
}

# If we're passed an experiment type we don't recognise, die.
else {
	die( "Don't know what to do with experiment of type \"", $atlasExperimentType, "\"\n" );
}


###############
# subroutines #
###############

# init
# 	Get commandline arguments using Getopt::Long.
# 		--atlasxml <Atlas contrasts XML file> (required)
# 		--irapconf <iRAP config file> (for RNA-seq)
sub init {
	
	use vars qw/ %opt /;
	use Getopt::Long;

	my ($atlasXML, $irapConfig);
	GetOptions(
		'atlasxml=s' => \$atlasXML,
		'irapconf=s' => \$irapConfig,
	);

	unless($atlasXML) {
		die "Cannot proceed without Atlas contrasts XML file. Please provide it with \"--atlasxml <filename>\"\nFor RNA-seq experiments please also provide \"--irapconf <iRAP config filename>\"\n";
	}
	
	return ($atlasXML, $irapConfig);
}


# run_microarray_differential_expression
# 	- For each array design, for each contrast, run differential expression
# 	analysis using limma script.
# 	- For each array design, combine results for all contrasts and write
# 	adjusted p-values, moderated t-statistics, and log2(fold-change)s to one
# 	file.
# 	- Create MvA plots for each contrast using MvA plotting script.
sub run_microarray_differential_expression {

	my ( $expAcc, $experimentConfig ) = @_;

	my $atlasProcessingDirectory = Cwd::cwd();

	# Ensure that XML config file exists in processing directory.
	my $configFile = File::Spec->catfile( $atlasProcessingDirectory, $expAcc . "-configuration.xml" );
	unless( -r $configFile ) {
		die( "Cannot find XML config file for $expAcc in $atlasProcessingDirectory\n" );
	}
	
	print( "Running differential expression analysis in R...\n" );

	# Run R script.
	my $R_limmaOutput = `$limmaScript $expAcc $atlasProcessingDirectory`;

	# Check R output for errors.
	if( $R_limmaOutput =~ /error/i ) {
		
		# Can't continue without results from limma so may as well quit.
		die("\nError during differential expression analysis, outout from R below.\n------------\n$R_limmaOutput\n------------\n");
	 
	} else {
		print "Differential expression analysis successful\n";
	}
			
	# Get the differential expression results into the hash of all analytics results.
	# Files are read out of the users ~/tmp directory where the R script wrote
	# them to.
	my $analyticsDEResults = read_limma_results( $expAcc );

	# Map contrast IDs to contrast names.
	my $contrastIDs2names = map_contrast_ids_to_names( $experimentConfig );
	my $contrastIDs2arrayDesigns = map_contrast_ids_to_arraydesigns( $experimentConfig );

	my $tempDir = File::Spec->catdir( $ENV{ "HOME" }, "tmp" );
	
	# Get the names of the MvA plot data files.
	my @plotDataFiles = glob( "$tempDir/$expAcc.g*_g*.plotdata.tsv" );

	foreach my $plotDataFile ( @plotDataFiles ) {

		( my $contrastID = basename( $plotDataFile ) ) =~ s/.*\.(g\d+_d\d+)\.plotdata.*/$1/;

		my $contrastName = $contrastIDs2names->{ $contrastID };

		# Filename for MvA plot.
		my $plotFile = $expAcc."_".$arrayDesignAccession."-".$contrastID."-mvaPlot.png";

		# TODO: get contrast names
		# Create MvA.
		make_mva_plot( "microarray", $plotFile, $plotDataFile, $contrastName, $mvaScript );

		`rm $plotDataFile`;
	}

	# Now we have results for all the contrasts in this analytics element. Write them to a file.
	write_results( $analyticsDEResults, $expAcc, $experimentConfig );
	
}


sub map_contrast_ids_to_names {

	my ( $experimentConfig ) = @_;

	my $allAnalytics = $experimentConfig->get_atlas_analytics;

	my $contrastIDs2arrayDesigns = {};

	foreach my $analytics ( @{ $allAnalytics } ) {

		my $contrasts = $analytics->get_atlas_contrasts;
		my $plaform = $analytics->get_platform;

		foreach my $contrast ( @{ $contrasts } ) {

			my $contrastID = $contrast->get_contrast_id;

			$contrastIDs2arrayDesigns->{ $contrastID } = $arrayDesign;
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


# join_assay_names
# 	- Join assay names using "<SEP>" ( or e.g. "<T1_SEP>" for technical
# 	replicates ).
sub join_assay_names {

	my ( $assays ) = @_;

	# First map assays to their technical replicate IDs (or
	# "no_technical_replicate_id" if they don't have one).
	my $techRepIDsToAssays = map_technical_replicate_ids_to_assays( $assays );

	# Array for assay names to join with "<SEP>"
	my $assayNamesToJoin = [];

	# Now go through the technical replicate IDs...
	foreach my $techRepID ( keys %{ $techRepIDsToAssays } ) {

		# For assays that aren't technical replicates...
		if( $techRepID eq "no_technical_replicate_id" ) {
			
			# Go through the assays
			foreach my $assay ( @{ $techRepIDsToAssays->{ $techRepID } } ) {
				
				# Get the assay name
				my $assayName = $assay->get_name;
				
				# Add it to the array of assay names to join.
				push @{ $assayNamesToJoin }, $assayName;
			}
		}
		# For assays that are technical replicates...
		else {

			# Create separator from technical replicate group ID, e.g.
			# "<T1_SEP>", "T2_SEP>", ...
			my $techRepSeparator = "<" . uc( $techRepID ) . "_SEP>";

			# Array for assay names for this technical replicate group.
			my $techRepAssays = [];

			# Get all the assay names for this technical replicate group.
			foreach my $assay ( @{ $techRepIDsToAssays->{ $techRepID } } ) {

				# Get the assay name.
				my $assayName = $assay->get_name;

				# Add it to the array of assay names for this technical replicate group.
				push @{ $techRepAssays }, $assayName;
			}

			# Now join the assay names with the separator.
			my $joinedTechRepAssayNames = join $techRepSeparator, @{ $techRepAssays };

			# Add this to the array of assay names to be joined
			push @{ $assayNamesToJoin }, $joinedTechRepAssayNames;
		}
	}

	# Now we've filled the array of assay names, with and without technical
	# replicates. Join them together with "<SEP>".
	my $joinedAssayNames = join "<SEP>", @{ $assayNamesToJoin };

	return $joinedAssayNames;
}


# read_limma_results
# 	- parses output from limma script and adds them to a hash.
sub read_limma_results {

	my ( $expAcc ) = @_;

	my $analyticsDEResults = {};

	my $tempDir = File::Spec->catdir( $ENV{ "HOME" }, "tmp" );

	my @limmaResultsFiles = glob( "$tempDir/$expAcc.g*_g*.analytics.tsv" );

	foreach my $limmaResultsFile ( @limmaResultsFiles ) {

		( my $contrastID = basename( $limmaResultsFile ) ) =~ s/.*\.(g\d+_d\d+)\.analytics.*/$1/;
	
		# Add results to hash.
		open(LIMMARES, "<$limmaResultsFile") or die("Can't open file $limmaResultsFile: $!\n");
		while(defined(my $line = <LIMMARES>)) {
			
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
		close(LIMMARES);
		# clean up
		`rm $limmaResultsFile`;
	}

	return $analyticsDEResults;
}


# strip_cy
# 	- remove ".Cy3" or ".Cy5" from ends of assay names.
sub strip_cy {
				
	my @assayNames = @_;
	my @assayNamesNoCy = ();
	foreach my $assayName (@assayNames) {

		$assayName =~ s/\.Cy\d$//;
		push @assayNamesNoCy, $assayName;
	}
	return @assayNamesNoCy;
}


# read_irap_conf
#  - read iRAP config file to get filename for DESeq results for each contrast.
#  - check experiment accession in iRAP config matches that of Atlas contrasts XML.
#  - add filename for DESeq results file for each contrast to $contrastHash.
sub read_irap_conf {

	my ( $irapConfig, $expAcc ) = @_;

	# Directory of config file is the same as the beginning of path to
	# DESeq files.
	my $irapDir = ( fileparse( $irapConfig ) )[1];

	# Get experiment accession, mapper, quantification method and DE method
	# from config file to create full path to DESeq results files.
	my ( $irapExptAcc, $mapper, $qMethod, $deMethod );
	
	open( CONF, $irapConfig ) or die( "Can't open $irapConfig: $!\n" );

	print "\nReading iRAP config: $irapConfig...\n";
	while( defined( my $line = <CONF> ) ) {

		chomp $line;

		if( $line =~ /^name=(.*)/ ) { 
			$irapExptAcc = $1; 
		
			# Check experiment accession from iRAP config matches that in filename of
			# Atlas contrasts file.
			unless( $irapExptAcc eq $expAcc ) {
				die "\niRAP experiment accession is: $irapExptAcc -- does not match Atlas experiment accession ($expAcc)\n";
			}
		}
		
		elsif( $line =~ /^mapper=(.*)/ ) { $mapper = $1; }
		elsif( $line =~ /^quant_method=(.*)/ ) { $qMethod = $1; }
		elsif( $line =~ /^de_method=(.*)/ ) { $deMethod = $1; }
	}
	close( CONF );

	# Build path to DESeq results directory.
	my $DESeqDirPath = File::Spec->catdir( $irapDir, $expAcc, $mapper, $qMethod, $deMethod );

	# Get all the DESeq results files.
	my @DESeqResultsFiles = glob( "$DESeqDirPath/*.genes_de.tsv" );

	# Empty hashref for DESeq results files.
	my $contrastIDsToDESeqFiles = {};

	# Go through the DESeq results files...
	foreach my $deseqFile ( @DESeqResultsFiles ) {
		
		# Get the contrast ID from the file name.
		( my $contrastID = basename( $deseqFile ) ) =~ s/\.genes_de\.tsv//;
		
		# Add the full path to the file to the hash under the contrast ID.
		$contrastIDsToDESeqFiles->{ $contrastID } = $deseqFile;
	}

	return ( $contrastIDsToDESeqFiles );
}


# get_rnaseq_results
# 	- For each contrast, parse the DESeq results file to get log2(fold-change)s
# 	and adjusted p-values.
# 	- Write all contrast results into one single file.
# 	- For each contrast, create MvA plots using MvA plotting script.
sub get_rnaseq_results {

	my ( $contrastIDsToDESeqFiles, $experimentConfig ) = @_;

	# Get the RNA-seq analytics section from the experiment config.
	my $allAnalytics = $experimentConfig->get_atlas_analytics;
	
	# Empty variable to put RNA-seq analytics object in.
	my $rnaSeqAnalytics;
	
	# Go through the analytics elements, get the rnaseq one.
	foreach my $analytics ( @{ $allAnalytics } ) {
		
		# Get the contrasts from the rnaseq analytics element.
		if( $analytics->get_platform eq "rnaseq" ) { $rnaSeqAnalytics = $analytics; }
	}	

	# Get the contrasts from the RNA-seq analytics object.
	my $rnaSeqContrasts = $rnaSeqAnalytics->get_atlas_contrasts;

	# Make a hash for the contrasts mapping IDs to names.
	my $contrastIDsToNames = {};
	foreach my $contrast ( @{ $rnaSeqContrasts } ) {

		my $contrastID = $contrast->get_contrast_id;
		my $contrastName = $contrast->get_contrast_name;

		$contrastIDsToNames->{ $contrastID } = $contrastName;
	}

	# Empty hash to store differential expression results before writing.
	my $rnaSeqDEResults = {};

	print "\nCollecting DESeq results...\n";
	
	# Go through the contrasts...
	foreach my $contrastID (keys %{ $contrastIDsToDESeqFiles } ) {

		print "Collecting results for contrast: ", $contrastIDsToNames->{ $contrastID }, "\n";
		
		# Get the path to the DESeq results file.
		my $deseqFile = $contrastIDsToDESeqFiles->{ $contrastID };

		# Need to get column indices of "baseMean", "log2FoldChange", and "padj" in DESeq results.
		my ($basemeanIdx, $logfcIdx, $adjpvalIdx);

		open(DESEQRES, $deseqFile) or die("Can't open $deseqFile: $!\n");

		# Temp file for MvA plot data.
		my $plotDataTempFile = "/tmp/plotData.$$.txt";
		
		open(PLOTDATA, ">$plotDataTempFile") or die("Can't open file to write temporary plot data to: $!\n");
		printf(PLOTDATA "ID\tavgExpr\tlogFC\tadjPval");

		print "Reading DESeq results from $deseqFile...";
		while(defined(my $line = <DESEQRES>)) {

			chomp $line;
			my @lineSplit = split "\t", $line;
			
			# On first line, get column indices
			if($line =~ /^id\t/) {

				$basemeanIdx = (grep { $lineSplit[$_] eq "baseMean" } 0..$#lineSplit)[0];
				$logfcIdx = (grep { $lineSplit[$_] eq "log2FoldChange" } 0..$#lineSplit)[0];
				$adjpvalIdx = (grep { $lineSplit[$_] eq "padj" } 0..$#lineSplit)[0];

				if($basemeanIdx == 0 || $logfcIdx == 0 || $adjpvalIdx == 0) {
					die("Couldn't get column indices for all required columns.\n");
				}
			}
			else {
				
				# Use indices found above to get values (gene ID is always the
				# very first element).
				my $geneID = $lineSplit[0];

				my $baseMean = $lineSplit[$basemeanIdx];
				my $logFC = $lineSplit[$logfcIdx];
				my $adjPval = $lineSplit[$adjpvalIdx];
				# Ensure we get numbers for the things that are supposed to
				# be numbers (weird prob with DESeq output).
				unless(looks_like_number($baseMean)) {
					die "\nDid not get numeric value for baseMean:\nGene ID: $geneID\nbaseMean: $baseMean\n";
				}
				# Allow "NA" adjusted p-values through.
				unless(looks_like_number($adjPval) || $adjPval eq "NA") {
					die "\nDid not get numeric value for adjusted p-value:\nGene ID: $geneID\nadjusted p-value: $adjPval\n";
				}
				# Allow "NA" log fold-changes through.
				unless(looks_like_number($logFC) || $logFC eq "NA") {
					die "\nDid not get numeric value for log2FoldChange:\nGene ID: $geneID\nlog2FoldChange: $logFC\n";
				}
				
				# Add to hash for file of all contrasts' results.
				$rnaSeqDEResults->{ $geneID }->{ $contrastID } = [ $adjPval, $logFC ];
				# Add to file for MvA plot.
				printf(PLOTDATA "\n$geneID\t$baseMean\t$logFC\t$adjPval");
			}
		}
		close(PLOTDATA);
		close(DESEQRES);
		print "done\n";

		# Filename for MvA plot
		my $plotFile = "$expAcc-".$contrastID."-mvaPlot.png";

		# Get the human-readable contrast name from the hash.
		my $atlasName = $contrastIDsToNames->{ $contrastID };

		# Create MvA
		make_mva_plot("rnaseq", $plotFile, $plotDataTempFile, $atlasName, $mvaScript);
	}
	
	# Write the results to a file.
	write_results( $rnaSeqDEResults, $experimentConfig->get_experiment_accession, $experimentConfig );
}


# make_mva_plot
# 	- Create MvA plot using MvA plot script.
sub make_mva_plot {

	$_ = shift for my ($platform, $plotFile, $plotDataTempFile, $contrastName, $mvaScript);

	print "Making MvA plot...";
	# Create MvA plot with MvA plot script
	my $R_mvaOutput = `$mvaScript $plotDataTempFile \"$contrastName\" $plotFile $platform 2>&1`;
	
	# Check for errors (this doesn't catch warnings).
	if($R_mvaOutput =~ /error/i) {

		# Report the error but don't worry about dying as we can live
		# without the MvA plot.
		print "\nError creating MvA plot, output from R below.\n------------\n$R_mvaOutput\n------------\n";

	} else {
		print "done.\n";
	}
	# clean up
	`rm $plotDataTempFile`;
	`rm Rplots.pdf`;
}


# write_results
#	- write all results for a single platform (e.g. "A-AFFY-35" or "rnaseq") to
#	a tab-delimited text file.
sub write_results {

	$_ = shift for my ( $analyticsDEResults, $experimentAccession, $experimentConfig);
	
	foreach my $analytics ( @{ $experimentConfig->get_atlas_analytics } ) {

		# Get the platform from the analytics element.
		my $platform = $analytics->get_platform;

		# Results file needs array design accession if it's microarray.
		my $resFile = ( $platform eq "rnaseq" ? $expAcc."-analytics.tsv.undecorated" : $expAcc."_".$platform."-analytics.tsv.undecorated" );

		# Get all the contrast IDs for this analytics element.
		my $analyticsContrastIDs = [];
		
		foreach my $contrast ( @{ $analytics->get_atlas_contrasts } ) {

			push @{ $analyticsContrastIDs }, $contrast->get_contrast_id;
		}
		
		print "\nWriting results to $resFile\n";

		open(RESFILE, ">$resFile") or die "Can't open $resFile: $!\n";

		# Write column headers
		if($platform eq "rnaseq") { printf(RESFILE "Gene ID"); }
		else { printf(RESFILE "Design Element"); }
		foreach my $contrastID ( @{ $analyticsContrastIDs }) {
			printf(RESFILE "\t$contrastID.p-value");
			unless($platform eq "rnaseq") { printf(RESFILE "\t$contrastID.t-statistic"); }
			printf(RESFILE "\t$contrastID.log2foldchange");
		}

		# Write statistics
		foreach my $id (keys %{ $analyticsDEResults }) {

			printf(RESFILE "\n$id");

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
				printf(RESFILE "\t$statsString");
			}
		}
		close( RESFILE );
	}
}
			

