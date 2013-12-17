#!/usr/bin/env perl

# POD documentation - main docs before the code
=pod

=head1 NAME

	Karyn Megy - 18-June-13
	kmegy@ebi.ac.uk
	diffAtlas_generateContrastsForExperiment.pl

=head1 SYNOPSIS

	Generate thecontrast file (XML) from a IDF/SDRF file (ArrayExpress format)

=head1 DESCRIPTION

	Generate the contrast file (XML) from a IDF/SDRF file (ArrayExpress format) and a config file. 

	The config file contains the list of reference terms and factor values to exclude from further analysis.
	It is formatted as follow:
		[REFERENCE] 
			wild type
			control
			no treatment 
			....
		[/REFERENCE]

		[FACTOR_VALUE_KILL]
			individual
			....
		[/FACTOR_VALUE_KILL]


=head1 OPTIONS

	See help 
	diffAtlas_generateContrastsForExperiment.pl -h

=head1 EXAMPLES

	diffAtlas_generateContrastsForExperiment.pl -exp experiment_name -conf config_file.txt -out path/to/output/XML/dir

	E.g. 
	   diffAtlas_generateContrastsForExperiment.pl -exp E-MTAB-1066  -conf analysis/differential/reference_assay_group_factor_values.txt
     	or 
       diffAtlas_generateContrastsForExperiment.pl -exp E-MTAB-1066  -conf analysis/differential/reference_assay_group_factor_values.txt -outdir /homes/kmegy/tmp 


=cut


use strict ;
use Getopt::Long qw(:config no_ignore_case) ;
use Magetab4Atlas ;


## Initialise global $, @ and %
my ($experiment, $conf, $referenceArg, $killFactorValue, $help, $debug) ; #arguments
my ($subDirectory, $commandLine) ; #values infered from command line
my $outdir = "." ; #default output directory
my %H_config ; #contains info. from the config file 
my %H_hybNameFactorValue ; #store association factor value / hybridization name
my @A_assayGroups ; #store assay groups 

my $errorCode = 0 ; #report error when testing the replicates, reference etc. 
my $errorMessage ; #error message when testing the replicates, reference etc.
my $flag ; #to parse config file
my $reference ; #reference Factor Value(s) to calculate D.E. against
my $noReferenceError ; # factor: candidate values - to output in log when reporting that no reference could be found     


## Get arguments
################
GetOptions( 
	'help|Help|h|H' => \$help,
	'exp=s'	 => \$experiment,
	'conf=s' => \$conf,
	'ref=s'  => \$referenceArg,
	'kill=s' => \$killFactorValue,
	'out=s'  => \$outdir,
	'debug=s'=> \$debug,
) ;

$commandLine = join(' ',@ARGV); 

if (!$experiment || !$conf) { print "[WARNING] Missing experiment (-exp $experiment) and/or configuration files (-conf $conf)\n" ; $help  = 1 ; }

if (!$outdir) { print "[WARNING] No output directory provided (-out $outdir). Default is current directory.\n" } ;

if ($help) { usage($commandLine) ; die ; }


## Directory with SDRF file
my $experimentDirectory = "/nfs/ma/home/arrayexpress/ae2_production/data/EXPERIMENT/" ;

## Experiment sub-directory
if ($experiment =~ /E-(\w+?)-\d+?/) { $subDirectory = $1 ; }
else { die "[ERROR] Experiment: $experiment not formatted as expected.\n" ; }

## IDF & SDRF files (input) and XML file (output)
my $idf = "$experimentDirectory/$subDirectory/$experiment/$experiment.idf.txt" ;
my $outfileXML = "$outdir/$experiment-configuration.xml.auto" ;


## Extract information from the config file
###########################################
open (CONF, $conf) || die "Can't open the configuration file $conf!\n" ;
while (my $line=<CONF>) {
	chomp $line ; 

	if ($line !~ /^#/) {

		#If reference value(s) or kill terms given in argument,
		#they take precedence over the ones from the config file.
		#Should be in double quotes and comma separated
		if ($referenceArg) {
				my @A_referenceArg = split(",", $referenceArg) ; 
				for my $refArg (@A_referenceArg) {
					#For each value, trim any start/end spaces (but not middle ones!)
					#and store
					$refArg =~ s/^\s+//g ; $refArg =~ s/\s+$//g ;
					$H_config{"REFERENCE"}{lc($refArg)} = 1 ; 
				}	
		}	
		if ($killFactorValue) { 	
			my @A_killFV = split(",", $killFactorValue) ;
			for my $factorValue (@A_killFV) {
				#For each value, trim any start/end spaces (but not middle ones!)
				#and store
				$factorValue =~ s/^\s+//g ; $factorValue =~ s/\s+$//g ;
				$H_config{"FACTOR_VALUE_KILL"}{lc($factorValue)} = 1 ; 
			}
		}

		#If no reference value(s) or kill terms  given in argument,
		#Take them from the config file.
		if ($line =~ /REFERENCE/ && !$referenceArg) { $flag = "REFERENCE" ; }
		elsif ($line =~ /FACTOR_VALUE_KILL/ && !$killFactorValue) { $flag = "FACTOR_VALUE_KILL" ; }
		elsif ($line =~ /\[\//) { $flag = "" ; }
		else { if (($flag ne "") && ($line ne "")) { $H_config{$flag}{lc($line)} = 1 ; } } #Use lc for case insensitive conparison
	}
}	
close CONF ;


## Collect FactorValues & ENA IDs
# from SDRF file. Use magetab2atlas module
#################################
# Using readMageta
if ($debug) { print "[DEBUG] Reading Magetab files\n" ; }
my ($factorType, $Href_efvs2runAccessions) = &readMagetab($idf) ;
my %H_eFactorValues2runIDs = %$Href_efvs2runAccessions ; #dereference the hash

## Factor Value parsing
## For each organism and each array design (usually: 1 of each only),
## Parse the factor values: 
#	- >= 3 replicates	      -- delete the Factor Value if not true 
#	- reference		          -- from the list of references generated from the config file, or passed as argument
#	- forbidden Factor Values -- from the kill list generated from the config file, or passed as argument; delete Factor Value if true

#Open the output XML file
open (XML, ">$outfileXML") || die ("Can't open output XML file $outfileXML\n") ;
my $configurationTag = 0 ;
$noReferenceError = "Candidate reference values for $factorType: ";
if ($debug) { print "[DEBUG] Parsing values collected in Magetab module\n" ; }
foreach my $species (keys %H_eFactorValues2runIDs) {
	if ($debug) { print "[DEBUG] Species is $species ($factorType)\n" ; }
	foreach my $array (keys %{$H_eFactorValues2runIDs{$species}}) {
		if ($debug) { print "[DEBUG]\tArray is $array\n" ; }

		foreach my $FV (keys %{$H_eFactorValues2runIDs{$species}{$array}}) {
			if ($debug) { print "[DEBUG] Testing '$FV' -- @{$H_eFactorValues2runIDs{$species}{$array}{$FV}}\n" ; }
			
			#Test for forbidden factor value (e.g. 'individual')
			if (exists $H_config{"FACTOR_VALUE_KILL"}{lc($FV)}) { delete $H_eFactorValues2runIDs{$species}{$array}{$FV} ; next ; } 		
			$noReferenceError .= " '$FV' ";

			#Test for reference - if already one, die loudly
			#(case insensitive: lc only)
			if (exists $H_config{"REFERENCE"}{lc($FV)}) { 
				if ($reference eq "") { $reference = $FV ; }
				else { die "[ERROR] $experiment - More than one reference: $reference and $FV!\n" } ;
			}

			#Test for replicates
			my $replicateCount = scalar @{$H_eFactorValues2runIDs{$species}{$array}{$FV}} ;
			if ($debug) { print "[DEBUG] Replicate number: $replicateCount\n" ; }
			if ($replicateCount < 3) { delete $H_eFactorValues2runIDs{$species}{$array}{$FV} ; }
		}	

		#Anything left afterwards
		print "[INFO] Checking Factor Values suitability for differential expression analysis\n" ;

		#Reference Factor Value ? 
		if (!defined $reference) { 
			$errorCode = 1 ;
			$errorMessage .= "No reference: $noReferenceError. " ;
		}

		#Any Factor Value left (>= 3 replicates)?
		#Need at least 2 of them!		
		if (keys %{$H_eFactorValues2runIDs{$species}{$array}} < 2) {
		#if (keys %H_hybNameFactorValue < 2) {
			$errorCode = 1 ;
			$errorMessage .= "Less than 2 values with at least 3 replicates. " ; 
		}	

		#If no error reported, then FactorValue test passed!
		if ($errorCode == 0) { 

			print "[INFO] Factor Value test passed!\n" ; 

			##Make groups 
			# Easier than making them on the fly when generating the XML
			# For simplicity, g1 will always be the reference and the rest assigned at random
			my $groupCounter = 1 ;
			$A_assayGroups[$groupCounter] = $reference ; #$A_assayGroups[0] empty, so that array position serves as group ID (g1, g2, g3 etc.) 

			foreach my $FV (keys %{$H_eFactorValues2runIDs{$species}{$array}}) {

				if ($FV ne $reference) {
					$groupCounter++ ;
					$A_assayGroups[$groupCounter] = $FV ;
				}		
			}

			print "[INFO] Print XML contrast file $outfileXML\n" ;
			##Format in XML
			#Beginning XML
			#If 1st one, print the 'configuration' tag
			if ($configurationTag == 0) { &XMLboundaries("start-conf") ; $configurationTag = 1 ; }
			&XMLboundaries("start-analytics") ; 

			#Array design section - only if experiment is an array
			if ($array ne "0") { &printArrayDesign($array) ; }	

			#Assay group section
			&printAssayGroup($species,$array) ;

			#Contrast section
			&printContrast($array) ;

			#End XML
			&XMLboundaries("end-analytics") ;

		} else {  #Cannot generate contrast file
			die "\n[ERROR] Contrast file cannot be generated for $experiment :: $species :: $array: $errorMessage\n" ; 
		}
	}		 	
}
if ($debug) { print "[DEBUG] Finish reading Magetab files\n" ; }

#If anything has been written in the XML contrast file
#Close the <configuration> tag
if ($configurationTag == 1) { &XMLboundaries("end-conf") ; }

#Close XML file
close XML ;
print "\n\n\n" ;



## Subroutine
#############
#Print usage for the program
sub usage {
	my ($command_line) = @_;
	
	print "Your command line was:\t".
		"$0 $command_line\n".
		"Compulsory parameters:\n".
		"\t-exp: experiment name. E.g. E-MTAB-123\n".
		"\t-conf: generic configuration file for that program. It should be in the same directory, and possible named reference_assay_group_factor_values.txt\n".
		"Optional parameters:\n".
		"\t-ref: list of possible reference terms to search for. In double quotes and comma separated if multiple. Take precedence over the config file.\n".
		"\t-kill: list of FactorValue terms to discard (to kill). In double quotes and comma separated if multiple. Take precedence over the config file.\n".
		"\t-outdir: output directory. Default is the current directory.\n" ;
}


#Print tabulations n times, in XML output file 
#'n' being given in argument
sub tabulationXML {
	print XML "\t" x $_[0] ;
}

# Print 'array_design' section
#Parameter is array name
sub printArrayDesign {
	my $subArray = $_[0] ;

	&tabulationXML(2) ; print XML "<array_design>\n" ;
	&tabulationXML(3) ; print XML "$subArray\n" ;
	&tabulationXML(2) ; print XML "</array_design>\n" ;
}

#Print 'assay_group' section
#Parameters are species and array names
sub printAssayGroup {

	my $subSpecies = $_[0] ;	
	my $subArray = $_[1] ;

	&tabulationXML(2) ; print XML "<assay_groups>\n" ;
	foreach my $i (1..$#A_assayGroups) { #there is nothing in $A_assayGroups[0]
		my $factVal = $A_assayGroups[$i] ;

		&tabulationXML(3) ; print XML "<assay_group id=\"g$i\">\n" ;
		foreach my $names (@{$H_eFactorValues2runIDs{$subSpecies}{$subArray}{$factVal}}) {
           	
			&tabulationXML(4) ;
			print XML "<assay>$names</assay>\n" ;
		}
		&tabulationXML(3) ; print XML "</assay_group>\n" ;
	}
	&tabulationXML(2) ; print XML "</assay_groups>\n" ;
}


#Print 'contrast' section
#Parameter is array name
sub printContrast {
	my $subArray = $_[0] ;

	$factorType = lc($factorType) ;

	&tabulationXML(2) ; print XML "<contrasts>\n" ;
	foreach my $i (2..$#A_assayGroups) { #starting at 2 because we have g1 (reference) vs. the rest
		
		&tabulationXML(3) ; print XML "<contrast id=\"g1_g".$i."\">\n" ;
		&tabulationXML(4) ; print XML "<name>$factorType:'$A_assayGroups[$i]' vs '$A_assayGroups[1]'" ;
		if ($subArray != 0) { print XML " on $subArray" ; } 
		print XML "</name>\n" ;
		&tabulationXML(4) ; print XML "<reference_assay_group>g1</reference_assay_group>\n" ;
		&tabulationXML(4) ; print XML "<test_assay_group>g".$i."</test_assay_group>\n" ;
		&tabulationXML(3) ; print XML "</contrast>\n" ;
	}	
	&tabulationXML(2) ; print XML "</contrasts>\n" ;
}	


#Print the beginning/end of the XML file
#Both are in the same subroutine to make it easier 
#to check that what's been open is being closed
sub XMLboundaries {
	my $location = $_[0] ;

	if ($location =~ /start-conf/) {
		print XML "<configuration>\n" ;
	}	
	
	if ($location =~ /start-analytics/) {
		&tabulationXML(1) ; print XML "<analytics>\n" ;
	}

	if ($location =~ /end-analytics/) {
		&tabulationXML(1) ; print XML "</analytics>\n" ;
	}	

	if ($location =~ /end-conf/) {
		print XML "</configuration>\n" ;
	}
}


#Read magetab file (SDRF / IDF files)
sub readMagetab {
    my $efvs2runAccessions = {} ;
	my $factorTypeString = "" ;

    if ($debug) { print "[DEBUG] readMagetab module\n" ; }

	# Create a Magetab4Atlas object. This reads the MAGETAB documents (via the
	# IDF filename provided) to get Atlas-relevant information for each assay in
	# the experiment.
	my $magetab4atlas = Magetab4Atlas->new( "idf_filename" => $idf );
	
	# Test if we find any assays.
	if(!$magetab4atlas->has_assays) { print "[DEBUG] No assay found!\n" ; } #die "No assay for this experiment!\n"; }
	 
	##Testing ...
	if ($debug) {print "[DEBUG] Experiment type: ".$magetab4atlas->get_experiment_type."\n"; }	
	
	if ($debug) {print "[DEBUG] Experiment type: ".$magetab4atlas->get_experiment_type."\n"; }

	my @A_magetabAssay = $magetab4atlas->get_assays ;
	if ($debug) {print "[DEBUG] Assays: $A_magetabAssay[0]\n"; }

	foreach my $assay4atlas (@{ $magetab4atlas->get_assays }) {

		#####
		## Issue: for some of the experiments, this doesn't return anything!
		# .... empty value!!
		#####
		if ($debug) { print "[DEBUG] Assays found!\n" ; }

		# Get the organism
		my $organism = $assay4atlas->get_organism() ; 
        print "\tOrganism:\n\t\t", $assay4atlas->get_organism, " ($organism) \n";

		# Get the assay name
		# Old SDRF: should be ENA_run (RNA-seq) or hybridisation name (microarray)
		# Newer SDRF: assay name field
		my $runAccession ;

		##if ($magetab4atlas->get_experiment_type =~ /array/) {
			$runAccession = $assay4atlas->get_name() ;
			print "Assay: ", $assay4atlas->get_name, " ($runAccession)\n";
		##}
		
		if ($magetab4atlas->get_experiment_type eq "RNA-seq") {

			if($assay4atlas->has_fastq_uri_set) {
				foreach my $ENArun (keys %{ $assay4atlas->get_fastq_uri_set }) {
					$runAccession = $ENArun ;
					print "\tENA run: $ENArun\n" ; 
				}
			}
		}	

		# Get the Factor type(s) and values for this assay
		print "\tFactors:\n";
		my $H_factors = $assay4atlas->get_factors;
		my $factorValueString = "" ; 
		$factorTypeString = "" ;
		foreach my $factorType (keys %{ $H_factors }) {

            $factorValueString .= $H_factors->{ $factorType }." " ;
			print "\t\t$factorType: ", $H_factors->{ $factorType }, "\n";

			#If >= 2 factor type, exclude the ones containing "*Time*"			
			unless ( keys %{$H_factors} >= 2 && $factorType =~ /Time/) { $factorTypeString .= $factorType ; }
		}

		$factorValueString =~ s/\s+?$// ; #remove trailing space
		print "\t\tfactorType string: $factorValueString\n";

		# Get array design for a microarray assay.
		my $arrayDesign = 0 ;
		if($magetab4atlas->get_experiment_type =~ /array/) {
			if($assay4atlas->has_array_design) {
				$arrayDesign = $assay4atlas->get_array_design() ;
				print "\tArray design:\n\t\t", $assay4atlas->get_array_design, " ($arrayDesign) \n" ;
			}
		}

		#Store
		if(exists($efvs2runAccessions->{ $organism }->{ $arrayDesign }->{ $factorValueString })) {
			push @{ $efvs2runAccessions->{ $organism }->{ $arrayDesign }->{ $factorValueString } }, $runAccession ;
		} else {
			$efvs2runAccessions->{ $organism }->{ $arrayDesign }->{ $factorValueString } = [ $runAccession ] ;
		}
	}

	#Return values
	#Return the factor value types (e.g. "genotype" - if multiple array: same F.V. types because it's the same SDRF file)
	# and the reference to a hash of mappings between factor values and run accessions
	return($factorTypeString, $efvs2runAccessions) ;

}
