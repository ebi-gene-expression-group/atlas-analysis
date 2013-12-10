#!/usr/bin/env perl
#
# Script for running microarray quality assessment. This script will read
# MAGE-TAB for an experiment and work out which data files belong to which
# array design, and are associated with which factor value(s). This info is
# written to a temp file which is read by an R script to do the quality
# assessment via the arrayQualityMetrics Bioconductor package.
#
# It will read the Atlas XML config defining assay groups and contrasts and
# only run QC for assays therein. It will remove assays that fail QC, and if
# the assay group(s) they belong to no longer have enough replicates as a
# result, it will remove those assay group(s) and the contrast(s) those assay
# group(s) are used in.
# If any assays fail QC, original XML config file will be renamed to:
# 	<expt accession>-configuration.xml.beforeQC
# New XML config will be written to:
# 	<expt accession>-configuration.xml

use strict;
use warnings;

# MAGE-TAB parsing.
use Magetab4Atlas;

# Helpful message
my $usage = "Usage:
	arrayQC.pl <experiment accession>
";

# Get accession from cmdline
my $exptAccession = shift;

# If nothing was provided, print message and die.
unless($exptAccession) { die $usage; }

# Make Atlas XML config filename. This script will run from the directory
# containing the XML file.
my $atlasXMLfile = "$exptAccession-configuration.xml";

# Die if the XML file doesn't exist.
unless(-e $atlasXMLfile) {
	die "[ERROR] Could not file $atlasXMLfile";
}

# Read XML config file and get hash of contrasts.
my ($H_contrastHash, $xmlExptType) = &readAtlasXML($atlasXMLfile);

# R script (should be in PATH).
my $qcRscript = "arrayQC.R";

# Path to directory with ArrayExpress/Atlas load directories underneath.
my $exptsLoadStem = "/ebi/microarray/home/arrayexpress/ae2_production/data/EXPERIMENT";

# miRBase mapped array designs -- we need to subset probes if we find one of these.
# Get an array of miRBase mapping files.
my @A_miRBaseFiles = glob("/ebi/microarray/home/atlas3-production/bioentity_properties/mirbase/*.A-*.tsv");

# Create a hash for easy checking.
my $H_miRBaseFileHash = {};
foreach my $miRBaseFile (@A_miRBaseFiles) {
	# Get the array design from the file name.
	(my $arrayDesign = $miRBaseFile) =~ s/.*(A-\w{4}-\d+)\.tsv/$1/;
	# Add the miRBase mapping file to the hash with the array design as key.
	$H_miRBaseFileHash->{ $arrayDesign } = $miRBaseFile;
}

# Get the pipeline (e.g. MEXP, MTAB, GEOD, ...) for this experiment.
(my $pipeline = $exptAccession) =~ s/E-(\w{4})-\d+/$1/;

# Experiment load directory and IDF filename.
my $loadDir = "$exptsLoadStem/$pipeline/$exptAccession";
my $idfFilename = "$loadDir/$exptAccession.idf.txt";

# Die if the IDF doesn't exist.
unless(-e $idfFilename) {
	die "[ERROR] Could not find IDF file $idfFilename\n";
}

# Read the MAGE-TAB.
my $magetab4atlas = Magetab4Atlas->new( "idf_filename" => $idfFilename );

# Next need to sort the raw data files by array design and within that
# by factor value. Use a hash like:
# 	$H->{ <array design 1> }->{ <factor value(s) 1> } = [ <file 1>, <file 2>, <file 3> ]
# 	                        ->{ <factor value(s) 2> } = [ <file 4>, <file 5>, <file 6> ]
# 	  ->{ <array design 2> }->{ <factor value(s) 3> } = [ <file 7>, <file 8>, <file 9> ]
# 	  ...
# Only consider assays that are in the XML config file (i.e. $H_contrastHash).
my ($H_arraysToFactorValuesToFiles, $experimentType) = &makeArraysToFactorValuesToFiles($magetab4atlas, $loadDir, $H_contrastHash);

# For each array design in the experiment: Write the factor value and filename
# information into a temporary text file for each array design. Run the R QC
# script and check output for errors. Remove assays that failed QC from the
# hash representing contrasts and assay groups, as well as assay group(s) and
# contrast(s) that are no longer eligible due to insufficient replication.
# Delete full paths to load directories from HTML report(s) produced by
# arrayQualityMetrics.
# Flag to set if any failed assays are found.
my $failed;
foreach my $arrayDesign (keys %{ $H_arraysToFactorValuesToFiles }) {
	
	# Write annotations (factor value(s), filenames, [labels]) to a temp file. This will be read by R and then deleted.
	my ($tempFile, $miRBaseFile) = &writeAnnotations($arrayDesign, $H_arraysToFactorValuesToFiles, $H_miRBaseFileHash, $experimentType);

	# Create directory name for this experiment/array design.
	my $reportDir = $exptAccession."_$arrayDesign"."_QM";
	# Run R script.
	my $qcRscriptOutput = `$qcRscript $tempFile $experimentType $exptAccession $arrayDesign $reportDir $miRBaseFile 2>&1`;

	# Check for errors in the R output.
	if($qcRscriptOutput =~ /error/i) {
		# Warn that QC had problems but continue with the next array design (if any).
		print "[QC] $exptAccession: Error during quality metrics calculation for array $arrayDesign, outout from R below.\n------------\n$qcRscriptOutput\n------------\n";
	}

	# Delete the no longer needed temp file.
	`rm $tempFile`;
	
	# Look for assays that failed QC and remove them from $H_contrastHash.
	($H_contrastHash, $failed) = &removeRejectedAssays($H_contrastHash, $qcRscriptOutput, $arrayDesign);

	# The HTML report contains the full path to the raw data files in th load
	# directory. This looks annoying and is not necessary, so remove it. The
	# path is in index.html and arrayQualityMetrics.js.
	&removeLoadDirFromReport($reportDir, $loadDir);
}

# Now check $H_contrastHash and see if it has any contrasts left. If not
# there's no point keeping this experiment at all. 
# First check each array design and remove the array design if there aren't any
# contrasts left.
foreach my $arrayDesign (keys %{ $H_contrastHash }) {
	# See if there are any keys left. If not, delete this array design.
	unless(keys %{ $H_contrastHash->{ $arrayDesign } }) {
		delete $H_contrastHash->{ $arrayDesign };
	}
}

# Now check whether whole experiment has any array designs left.
unless(keys %{ $H_contrastHash }) {
	# Log that there aren't any contrasts left to STDOUT.
	print "[QC] $exptAccession no longer has any eligible contrasts.\n";
	
	# Rename XML config file.
	print "\nRenaming $atlasXMLfile to $atlasXMLfile.beforeQC\n";
	`mv $atlasXMLfile $atlasXMLfile.beforeQC`;
	
	# Quit now because there's no point writing a new XML file. Exit 1 so that
	# we know there was an error.
	exit 1;
}

# Rewrite XML config file without failing assays and contrasts without enough
# replicates as a result. Use assay group IDs from contrast (assay group pair)
# IDs to write back to XML.
if($failed) {
	&writeNewXML($atlasXMLfile, $H_contrastHash);
}
# end
#####


###############
# Subroutines #
###############

# readAtlasXML
#	TODO: This is copied from diffAtlas_DE.pl -- should be in a module (ticket FGC-30).
# 	- read Atlas XML contrast definitions file.
#	- This file has an <analytics> section for each platform. For microarray
#	each analytics section has an array design accession in <array_design>
#	tags. We need to know which contrasts belong to which array design because
#	we have to select the appropriate normalized expressions file when we come
#	to do the DE analysis with limma.
#	- RNA-seq analytics section doesn't have an <array_design> section. All
#	RNA-seq data goes into the same <analytics> section. If there are multiple
#	reference genomes used this is taken care of in the iRAP config file but
#	this script doesn't care about that, because it is just retrieving
#	pre-computed results for each contrast found in the XML and iRAP config.
sub readAtlasXML {
	my ($atlasXML) = @_;
	
	# empty hash for contrast info.
	my $H_contrastHash = {};
	
	# load XML::Simple with strict mode -- gives helpful error messages.
	use XML::Simple qw(:strict);
	
	print "\nReading Atlas XML config from $atlasXML...\n";
	# xml contains assay group and contrast definitions parsed from XML file, for each <analytics> section.
	my $xml = XMLin($atlasXML, ForceArray => ['analytics', 'assay', 'contrast', 'assay_group'], KeyAttr => { contrast => 'id', assay_group => 'id' });
	
	my $xmlExptType = $xml->{ "experimentType" };
	
	# For each platform (array design or RNA-seq section)
	foreach my $ana (@{ $xml->{ "analytics" } }) { 
		my $arrayDesign = $ana->{ "array_design" };

		if(defined($arrayDesign)) { $arrayDesign =~ s/\s//g; }
		else { $arrayDesign = "rnaseq"; }

		my $thisContrast = $ana->{ "contrasts" }->{ "contrast" };
		foreach my $agPair (keys %{ $thisContrast }) {
			$H_contrastHash->{ $arrayDesign }->{ $agPair }->{ "atlasName" } = $thisContrast->{ $agPair }->{ "name" };
			unless($arrayDesign eq "rnaseq") {
				my $refAGnum = $thisContrast->{ $agPair }->{ "reference_assay_group" };
				my $testAGnum = $thisContrast->{ $agPair }->{ "test_assay_group" };

				# Have to create a new array for the assay accessions of each assay
				# group by de-referencing the arrayref in $ana.
				# Then pass the reference to this to $H_contrastHash.
				# This is because if the same assay group is used more than once,
				# passing the arrayref in $ana straight to $H_contrastHash does not
				# work. The first time it does, but the second time the assay group
				# is passed you just get e.g. (from Data::Dumper): 
				# 	'reference' => $VAR1->{'A-AFFY-35'}{'g2_g3'}{'reference'}
				# instead of:
				#	'reference' => [
				#                   'WT3',
				#                   'WT1',
				#                   'WT2'
				#                  ]
				# The second thing is what we actually wanted.
				# The following code achieves that.
				my @refAGarray = @{ $ana->{ "assay_groups" }->{ "assay_group" }->{ $refAGnum }->{ "assay" } };	
				$H_contrastHash->{ $arrayDesign }->{ $agPair }->{ "reference" } = \@refAGarray;

				my @testAGarray = @{ $ana->{ "assay_groups" }->{ "assay_group" }->{ $testAGnum }->{ "assay" }};
				$H_contrastHash->{ $arrayDesign }->{ $agPair }->{ "test" } = \@testAGarray;
			}
		}
	}

	# Log what we've found
	foreach my $arrayDesign (keys %{ $H_contrastHash }) {
		my $numContrasts = keys %{ $H_contrastHash->{ $arrayDesign }};
		print "$exptAccession: $numContrasts contrast";
		unless($numContrasts == 1) { print "s";} 
		print " found for $exptAccession - $arrayDesign:\n";
		
		foreach my $agPair (keys %{ $H_contrastHash->{ $arrayDesign } }) {

			print "\t", $H_contrastHash->{ $arrayDesign}->{ $agPair }->{ "atlasName" }, "\n";
		}
	}
	print "\n";
	return ($H_contrastHash, $xmlExptType);
}


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
# 	- $H_contrastHash : hash representation of all contrasts from Atlas XML config file.
sub makeArraysToFactorValuesToFiles {
	# Magetab4Atlas object and path to load directory.
	my ($magetab4atlas, $loadDir, $H_contrastHash) = @_;

	# Experiment type from Magetab4Atlas will be either "one-colour array" or
	# "two-colour array". Die if it's something else.
	my $experimentType = $magetab4atlas->get_experiment_type;
	unless($experimentType eq "one-colour array" || $experimentType eq "two-colour array") {
		die "This doesn't look like a microarray experiment. Experiment type found is: $experimentType\n";
	}

	# Get unique assay names from $H_contrastHash.
	my $H_xmlAssayNames = {};
	foreach my $arrayDesign (keys %{ $H_contrastHash }) {
		foreach my $assayGroupPair (keys %{ $H_contrastHash->{ $arrayDesign } }) {
			foreach my $testAssay (@{ $H_contrastHash->{ $arrayDesign }->{ $assayGroupPair }->{ "test" }}) {
				$H_xmlAssayNames->{ $testAssay } = 1;
			}
			foreach my $refAssay (@{ $H_contrastHash->{ $arrayDesign }->{ $assayGroupPair }->{ "reference" }}) {
				$H_xmlAssayNames->{ $refAssay } = 1;
			}
		}
	}
			
	# Ref to empty hash to fill.
	my $H_arraysToFactorValuesToFiles = {};
	# Go through the assays...
	foreach my $assay4atlas (@{ $magetab4atlas->get_assays }) {
		# Get the assay name
		my $assayName = $assay4atlas->get_name;
		# Skip this one if it's not in the XML config file
		unless(exists($H_xmlAssayNames->{ $assayName })) { next; }

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


# writeAnnotations
# 	- Writes the raw data filenames, factor values and assay names [and labels
# 	for 2-colour] to a temporary file.
# ARGUMENTS:
# 	- $arrayDesign : ArrayExpress array design accession
# 	- $H_arraysToFactorValuesToFiles : reference to hash of array designs, factor values, assay names and filenames.
# 	- $H_miRBaseFileHash : reference to hash of array designs that are for microRNA.
# 	- $experimentType : affy, agil1, or agil2.
sub writeAnnotations {
	my ($arrayDesign, $H_arraysToFactorValuesToFiles, $H_miRBaseFileHash, $experimentType) = @_;

	# Check if the array design has a miRBase mapping file
	# Flag
	my $miRBaseFile = 0;
	if(exists($H_miRBaseFileHash->{ $arrayDesign })) {
		$miRBaseFile = $H_miRBaseFileHash->{ $arrayDesign };
	}
	
	# Name for temp file to write annotations to, in /tmp with process ID ($$).
	my $tempFile = "/tmp/$exptAccession"."_$arrayDesign.$$.tsv";
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

	return ($tempFile, $miRBaseFile);
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
		# A filename for a temp file to write to, in /tmp with process ID.
		my $temp = "/tmp/$$.temp";
		
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


# &removeRejectedAssays
#	- Find names of assays that failed QC in the R script output and remove
#	them from $H_contrastHash. Also remove contrasts containing them if the
#	contrast no longer has enough replicates without failed assays.
# ARGUMENTS:
# 	- $H_contrastHash : reference to hash representation of Atlas XML config detailing all contrasts in experiment
# 	- $qcRscriptOutput : variable containing all output from R script (STDOUT, STDERR)
# 	- $arrayDesign : ArrayExpress array design accession
sub removeRejectedAssays {
	my ($H_contrastHash, $qcRscriptOutput, $arrayDesign) = @_;
	
	# Flag to set if we see any failed assays (hence XML needs re-write).
	my $failed = 0;

	# Check for rejected assays. If there are some, check them in the contrast hash.
	if($qcRscriptOutput =~ /REJECTED ASSAYS:\t(.*)\n/) {
		# Get the string of rejected assay names separated by tabs if more than one.
		my @A_rejectedAssays = split "\t", $1;

		# Want to check whether assay groups containing the rejected assay(s)
		# still have enough assays without rejected ones.
		foreach my $rejected (@A_rejectedAssays) { 
			# Log that this assay failed.
			print "[QC] $exptAccession: Assay \"$rejected\" failed QC and will be removed from XML config.\n";
			
			# Set flag
			$failed = 1;

			# Now see which contrasts it was found in and whether this means we
			# have to reject the contrast. If this contrast has to be rejected,
			# delete it from $H_contrastHash -- that way we can look at the end
			# and see if any contrasts are left in the experiment.
			foreach my $assayGroupPair (keys %{ $H_contrastHash->{ $arrayDesign } }) {
				my @A_testAssays = @{ $H_contrastHash->{ $arrayDesign }->{ $assayGroupPair }->{ "test" } };
				my @A_refAssays = @{ $H_contrastHash->{ $arrayDesign }->{ $assayGroupPair }->{ "reference" } };

				if(grep $_ eq $rejected, @A_testAssays) {
					# Get the contrast name.
					my $contrastName = $H_contrastHash->{ $arrayDesign }->{ $assayGroupPair }->{ "atlasName" };
					
					# Log the contrast name the assay was found in.
					print "[QC] $exptAccession: Assay \"$rejected\" found in test assay group in contrast \"$contrastName\".\n";

					# Make a new array without the rejected assay
					my @A_newTestAssays = ();
					foreach my $assayName (@A_testAssays) {
						unless($assayName eq $rejected) {
							push @A_newTestAssays, $assayName;
						}
					}
				 
					# Replace the array of test assays in $H_contrastHash with the new @A_newTestAssays
					$H_contrastHash->{ $arrayDesign }->{ $assayGroupPair }->{ "test" } = \@A_newTestAssays;

					# If there's only three assays in the assay group
					# containing the rejected one, removing it would mean too
					# few replicates. So delete this contrast from
					# $H_contrastHash.
					if(@A_newTestAssays < 3) {
						# Log that
						print "[QC] $exptAccession: Contrast \"$contrastName\" is no longer eligible: the test assay group no longer has enough replicates.\n";
						print "[QC] Removing contrast \"$contrastName\" from XML config.\n";
						
						# Remove this contrast from $H_contrastHash.
						delete $H_contrastHash->{ $arrayDesign }->{ $assayGroupPair };
					}
					else {
						print "Contrast \"$contrastName\" is still eligible.\n\n";
					}
				} 
				elsif(grep $_ eq $rejected, @A_refAssays) {
					# Get the contrast name.
					my $contrastName = $H_contrastHash->{ $arrayDesign }->{ $assayGroupPair }->{ "atlasName" };
				
					# Log the contrast name the assay was found in.
					print "[QC] $exptAccession: Assay \"$rejected\" found in reference assay group in contrast \"$contrastName\"\n";

					# Make a new array without the rejected assay
					my @A_newRefAssays = ();
					foreach my $assayName (@A_refAssays) {
						unless($assayName eq $rejected) {
							push @A_newRefAssays, $assayName;
						}
					}

					# Replace the array of reference assays in $H_contrastHash with the new @A_refAssays
					$H_contrastHash->{ $arrayDesign }->{ $assayGroupPair }->{ "reference" } = \@A_newRefAssays;
					
					# If there's only three assays in the assay group
					# containing the rejected one, removing it would mean too
					# few replicates. So delete this contrast from
					# $H_contrastHash.
					if(@A_newRefAssays < 3) {
						# Log that
						print "[QC] $exptAccession: Contrast \"$contrastName\" is no longer eligible: the reference assay group no longer has enough replicates\n\n";
						print "[QC] Removing contrast \"$contrastName\" from XML config.\n";
						
						# Remove this contrast from $H_contrastHash.
						delete $H_contrastHash->{ $arrayDesign }->{ $assayGroupPair };
					}
					else {
						print "Contrast \"$contrastName\" is still eligible.\n\n";
					}
				}
			}
		}
	}
	return($H_contrastHash, $failed);
}


# writeNewXML
# 	- Re-writes Atlas XML configuration file without assays that failed QC,
# 	using $H_contrastHash now it does not contain those assays.
# ARGUMENTS:
# 	- $atlasXMLfile : filename of original Atlas XML configuration.
# 	- $H_contrastHash : reference to hash representation of all contrasts and
# 	their assay groups that passed QC.
sub writeNewXML {
	my ($atlasXMLfile, $H_contrastHash) = @_;

	# Modules for writing XML.
	use XML::Writer;
	use IO::File;
	
	print "Writing new XML config.\n";

	# grep to find contrast IDs in original XML file.
	my $contrastTagLines = `grep \"<contrast id=\" $atlasXMLfile`;

	# Split the tags string on newlines.
	my @A_contrastTags = split "\n", $contrastTagLines;

	# Reference to empty array to store contrast IDs in the order they were in the
	# original XML file.
	my $A_contrastIDorder = [];

	foreach my $contrastTag (@A_contrastTags) {
		# Strip everything but the contrast ID.
		$contrastTag =~ s/.*id="(.*)">/$1/;
		# Add it to the array of contrast IDs to remember contrast order.
		push @{ $A_contrastIDorder }, $contrastTag;
	}
	
	# Now rename old XML config file.
	print "\nRenaming $atlasXMLfile to $atlasXMLfile.beforeQC\n";
	`mv $atlasXMLfile $atlasXMLfile.beforeQC`;

	# Open file to write to.
	my $newXML = IO::File->new(">$atlasXMLfile");
	
	# New XML writer that does newlines and nice indentation.
	my $writer = XML::Writer->new(OUTPUT => $newXML, DATA_MODE => 1, DATA_INDENT => 4);

	# Begin configuration XML, add experiment type.
	$writer->startTag("configuration", "experimentType" => $xmlExptType);

	foreach my $arrayDesign (keys %{ $H_contrastHash }) {
		# Start analytics element for this array design.
		$writer->startTag("analytics");

		# Add array design.
		$writer->dataElement("array_design", $arrayDesign);
		
		# Start assay_groups element for this array design.
		$writer->startTag("assay_groups");

		# Get assay groups for this array design.
		my $H_assayGroups = &getAssayGroups($H_contrastHash->{ $arrayDesign });

		# Go through the hash an add the assay groups to the XML.
		foreach my $assayGroupID (sort keys %{ $H_assayGroups }) {
			# Start the assay_group element for this assay group and add the ID.
			$writer->startTag("assay_group", "id" => $assayGroupID);
			
			# Go through the assays...
			foreach my $assay (@{ $H_assayGroups->{ $assayGroupID } }) {
				# Add assay elements to XML.
				$writer->dataElement("assay", $assay);
			}

			# Close assay_group element.
			$writer->endTag("assay_group");
		}

		# Close assay_groups element.
		$writer->endTag("assay_groups");

		# Start contrasts element for this array design.
		$writer->startTag("contrasts");

		# Add contrasts for this array design.
		# Go through the $A_contrastIDorder array to preserve ordering of contrasts
		# from original XML file.
		foreach my $assayGroupPair (@{ $A_contrastIDorder }) {
			# Check the contrast is stil in $H_contrastHash
			if(exists($H_contrastHash->{ $arrayDesign }->{ $assayGroupPair })) {
				# Start contrast element with ID.
				$writer->startTag("contrast", "id" => $assayGroupPair);
				
				# Add the contrast name element.
				$writer->dataElement("name", $H_contrastHash->{ $arrayDesign }->{ $assayGroupPair }->{ "atlasName" });

				# Add the reference and test assay group ID elements.
				my ($refGroupID, $testGroupID) = split "_", $assayGroupPair;
				$writer->dataElement("reference_assay_group", $refGroupID);
				$writer->dataElement("test_assay_group", $testGroupID);

				# Close contrast element.
				$writer->endTag("contrast")
			}
		}
		
		# Close contrasts element.
		$writer->endTag("contrasts");

		# Close analytics element.
		$writer->endTag("analytics");
	}

	# End configuration element.
	$writer->endTag("configuration");
	$writer->end;
	# Close file.
	$newXML->close;
}


# getAssayGroups
# 	- Takes assay groups and their IDs from a hash of contrasts for an array
# 	design, and returns a hash with assay group IDs as keys and references to
# 	arrays of assay groups as values.
sub getAssayGroups {
	my ($H_arrayDesignContrasts) = @_;

	# Reference to empty hash for assay groups.
	my $H_assayGroups = {};

	# Go through the contrasts for this array design...
	foreach my $assayGroupPair (keys %{ $H_arrayDesignContrasts }) {
		# Get the test and reference assay groups
		my $A_testAssays = $H_arrayDesignContrasts->{ $assayGroupPair }->{ "test" };
		my $A_refAssays = $H_arrayDesignContrasts->{ $assayGroupPair }->{ "reference" };

		# Get the reference and test assay group IDs from the contrast ID ($assayGroupPair).
		my ($refGroupID, $testGroupID) = split "_", $assayGroupPair;
		
		# If they're not already in the hash, add the assay groups with their IDs.
		unless(exists($H_assayGroups->{ $testGroupID })) {
			$H_assayGroups->{ $testGroupID } = $A_testAssays;
		}
		unless(exists($H_assayGroups->{ $refGroupID })) {
			$H_assayGroups->{ $refGroupID } = $A_refAssays;
		}
	}
	return $H_assayGroups;
}
