#!/usr/bin/env perl
#
# Script for running microarray quality assessment. This script will read
# MAGE-TAB for an experiment and work out which data files belong to which
# array design, and are associated with which factor value(s). This info is
# written to a temp file which is read by an R script to do the quality
# assessment via the arrayQualityMetrics Bioconductor package.

use strict;
use warnings;

# MAGE-TAB parsing.
use Magetab4Atlas;

my $exptAccession = shift;

my $qcRscript = "arrayQC.R";

# Path to directory with ArrayExpress/Atlas load directories underneath.
my $exptsLoadStem = "/nfs/ma/home/arrayexpress/ae2_production/data/EXPERIMENT";

# miRBase mapped array designs -- we need to subset probes if we find one of these.
my $miRBaseFiles = `ls /nfs/ma/home/atlas3-production/bioentity_properties/mirbase/*.A-*.tsv`;
# Put them in an array.
my @A_miRBaseFiles = split "\n", $miRBaseFiles;

# Get the pipeline (e.g. MEXP, MTAB, GEOD, ...) for this experiment.
(my $pipeline = $exptAccession) =~ s/E-(\w{4})-\d+/$1/;

# Experiment load directory and IDF filename.
my $loadDir = "$exptsLoadStem/$pipeline/$exptAccession";
my $idfFilename = "$loadDir/$exptAccession.idf.txt";

# Read the MAGE-TAB.
my $magetab4atlas = Magetab4Atlas->new( "idf_filename" => $idfFilename );

# Next need to sort the raw data files by array design and within that
# by factor value. Use a hash like:
# 	$H->{ <array design 1> }->{ <factor value(s) 1> } = [ <file 1>, <file 2>, <file 3> ]
# 	                        ->{ <factor value(s) 2> } = [ <file 4>, <file 5>, <file 6> ]
# 	  ->{ <array design 2> }->{ <factor value(s) 3> } = [ <file 7>, <file 8>, <file 9> ]
# 	  ...
my ($H_arraysToFactorValuesToFiles, $experimentType) = &makeArraysToFactorValuesToFiles($magetab4atlas, $loadDir);

# Write the factor value and filename information into a temporary text file
# for each array design. This will be read by R and then deleted.
foreach my $arrayDesign (keys %{ $H_arraysToFactorValuesToFiles }) {
	# Check if the array design has a miRBase mapping file
	# Flag
	my $miRBaseFile = 0;
	# Go through the files and look for a match.
	# This bit could end up taking a while if there are a lot of files?
	foreach my $file (@A_miRBaseFiles) {
		# If the array design matches log and remember the file.
		if($file =~ /$arrayDesign/) {
			print "miRNA array detected\n";
			$miRBaseFile = $file;
		}
	}
	
	# Name for temp file to write annotations to.
	my $tempFile = ".$exptAccession"."_$arrayDesign.tsv";
	# Create annotations temp file.
	open(my $tmpFH, '>', $tempFile) or die "Can't create file \"$tempFile\": $!\n";
	
	# Write factor values and corresponding raw data filename.
	# If this is two-colour data we want to include the label info in the file.
	# Create a file with headings like:
	# 	AssayName	Cy3	Cy5	FileName
	if($experimentType eq "agil2") {
		# Ref to empty hash to remember two-colour annotations and filenames.
		my $H_twoColourAnnotations = {};
		# Go through the factor values for this array design...
		foreach my $factorValue (keys %{ $H_arraysToFactorValuesToFiles->{ $arrayDesign } }) {
			# Go through the assays for this factor value...
			foreach my $assayName (keys %{ $H_arraysToFactorValuesToFiles->{ $arrayDesign }->{ $factorValue } }) {
				# Get the label (Cy3 or Cy5)
				(my $label = $assayName) =~ s/.*\.(Cy\d)$/$1/;
				# Make a version of the assay name without the label.
				(my $assayNameNoLabel = $assayName) =~ s/\.Cy\d$//;
				
				# Add the factor value to the hash for this assay for this label.
				$H_twoColourAnnotations->{ $assayNameNoLabel }->{ $label } = $factorValue;
				# Add the file name for this assay as well.
				$H_twoColourAnnotations->{ $assayNameNoLabel }->{ "filename" } = $H_arraysToFactorValuesToFiles->{ $arrayDesign }->{ $factorValue }->{ $assayName };
			}
		}

		# Write header.
		print $tmpFH "AssayName\tCy3\tCy5\tFileName";
		# Write the annotations.
		foreach my $assayNameNoLabel (keys %{ $H_twoColourAnnotations }) {
			print $tmpFH "\n$assayNameNoLabel\t";
			print $tmpFH $H_twoColourAnnotations->{ $assayNameNoLabel }->{ "Cy3" }, "\t";
			print $tmpFH $H_twoColourAnnotations->{ $assayNameNoLabel }->{ "Cy5" }, "\t";
			print $tmpFH $H_twoColourAnnotations->{ $assayNameNoLabel }->{ "filename" };
		}
	}
	# For 1-colour arrays, don't need the Cy3/Cy5 info.
	else {
		# Write header.
		print $tmpFH "AssayName\tFactorValue\tFileName";
		# Go through the factor values for this array...
		foreach my $factorValue (keys %{ $H_arraysToFactorValuesToFiles->{ $arrayDesign } }) {
			# Go through the assays for this factor value...
			foreach my $assayName (keys %{ $H_arraysToFactorValuesToFiles->{ $arrayDesign }->{ $factorValue } }) {
				# Write annotations.
				print $tmpFH "\n$assayName\t$factorValue\t".$H_arraysToFactorValuesToFiles->{ $arrayDesign }->{ $factorValue }->{ $assayName };
			}
		}
	}
	close $tmpFH;

	# Create directory name for this experiment/array design.
	my $reportDir = $exptAccession."_$arrayDesign"."_QM";
	# Run R script.
	my $qcRscriptOutput = `$qcRscript $tempFile $experimentType $exptAccession $arrayDesign $reportDir $miRBaseFile 2>&1`;

	# Check for errors in the R output.
	if($qcRscriptOutput =~ /error/i) {
		# Warn that QC had problems but continue with the next array design (if any).
		print "[WARNING] Error during quality metrics calculation for array $arrayDesign, outout from R below.\n------------\n$qcRscriptOutput\n------------\n";
	}
	
	# Delete the no longer needed temp file.
	`rm $tempFile`;

	# The HTML report contains the full path to the raw data files in th load
	# directory. This looks annoying and is not necessary, so remove it. The
	# path is in index.html and arrayQualityMetrics.js.
	&removeLoadDirFromReport($reportDir, $loadDir);
}


# Subroutines

# &makeArraysToFactoValuesToFiles
#  - Creates a hash sorting out the data files by array design and then factor value.
#  - E.g.:
# 	$H->{ <array design 1> }->{ <factor value(s) 1> }->{ <assay name 1> } = <file 1>
# 	                        ->{ <factor value(s) 2> }->{ <assay name 2> } = <file 2>
# 	  ->{ <array design 2> }->{ <factor value(s) 3> }->{ <assay name 3> } = <file 3>
# 	  ...
# Arguments:
# 	- $magetab4atlas : a Magetab4Atlas object
# 	- $loadDir : path to load directory containing raw data files.
sub makeArraysToFactorValuesToFiles {
	# Magetab4Atlas object and path to load directory.
	my ($magetab4atlas, $loadDir) = @_;
	
	# Experiment type from Magetab4Atlas will be either "one-colour array" or
	# "two-colour array". Die if it's something else.
	my $experimentType = $magetab4atlas->get_experiment_type;
	unless($experimentType eq "one-colour array" || $experimentType eq "two-colour array") {
		die "This doesn't look like a microarray experiment. Experiment type found is: $experimentType\n";
	}
	
	# Ref to empty hash to fill.
	my $H_arraysToFactorValuesToFiles = {};
	# Go through the assays...
	foreach my $assay4atlas (@{ $magetab4atlas->get_assays }) {
		# Get the assay name
		my $assayName = $assay4atlas->get_name;
		# Get the array design
		my $arrayDesign = $assay4atlas->get_array_design;
		# Get the factor(s) and value(s)
		my $H_factors = $assay4atlas->get_factors;
		# Get the raw data filename
		my $arrayDataFile = $loadDir."/".$assay4atlas->get_array_data_file;

		# For 1-colour array data, need to tell the R script whether this is
		# Affymetrix or Agilent (or other -- for now we only handle Affy and Agil
		# data). Work this out based on the file extension for now -- another
		# possibility might be to parse the ADF, but there is no standard way
		# to record manufacturer there so that might not be ideal.
		if($experimentType eq "one-colour array") {
			if($arrayDataFile =~ /\.cel$/i) { $experimentType = "affy"; }
			else { $experimentType = "agil1"; }
		} elsif($experimentType eq "two-colour array") {
			$experimentType = "agil2";
		}

		# Push all factor values onto an array.
		my @A_factorValues = ();
		foreach my $factor (keys %{ $H_factors }) {
			push @A_factorValues, $H_factors->{ $factor };
		}

		# Stick the factor values together if there's more than one. If there's
		# only one this just returns the factor value by itself.
		my $factorValue = join ", ", @A_factorValues;
		
		$H_arraysToFactorValuesToFiles->{ $arrayDesign }->{ $factorValue }->{ $assayName } = $arrayDataFile;
	}
	return ($H_arraysToFactorValuesToFiles, $experimentType);
}


# &removeLoadDirFromReport
# 	- Remove the full path to the load directory from index.html and arrayQualityMetrics.js.
# Arguments:
# 	- $reportDir : the directory containing the output arrayQualityMetrics.
# 	- $loadDir : the load directory path to remove.
sub removeLoadDirFromReport {
	my ($reportDir, $loadDir) = @_;
	
	# We need to fix the HTML report and the JavaScript file.
	my @A_filesToFix = ("$reportDir/index.html", "$reportDir/arrayQualityMetrics.js");

	foreach my $original (@A_filesToFix) {
		# A filename for a temp file to write to.
		my $temp = "$reportDir/.temp";
		
		# Open the report for reading.
		open(my $originalFH, "<", $original) or die "Can't open $original : $!\n";
		# Open the temp file for writing.
		open(my $tempFH, ">", $temp) or die "Can't open $temp : $!\n";
		
		# Go through the report line by line...
		while(defined(my $line = <$originalFH>)) {
			# Remove the load directory path if it's there.
			$line =~ s/$loadDir\/*//g;
			# Write the line to the temp file.
			print $tempFH $line;
		}
		# Close them.
		close($originalFH);
		close($tempFH);

		# Overwrite the original report with the fixed one.
		`mv $temp $original`;
	}
}

