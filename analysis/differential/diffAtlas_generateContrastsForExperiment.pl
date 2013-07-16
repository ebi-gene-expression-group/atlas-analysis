#!/usr/bin/env perl

# POD documentation - main docs before the code
=pod

=head1 NAME

  Karyn Megy - 18-June-13
  kmegy@ebi.ac.uk
  diffAtlas_generateContrastsForExperiment.pl

=head1 SYNOPSIS

  Generate a contrast file (XML) from SDRF file (ArrayExpress format)

=head1 DESCRIPTION

  Generate a contrast file (XML) from a SDRF file (ArrayExpress format) and a config file. 

  The config file contains the list of reference terms and factor values to exclude from further analysis.

  Note: 
	For now the SDRF parsing is with a gawk command
	Later: use Maria's code  


=head1 OPTIONS

  none

=head1 EXAMPLES

  diffAtlas_generateContrastsForExperiment.pl -exp experiment_name -conf config_file.txt -outdir path/to/output/XML/file

  E.g. 
	diffAtlas_generateContrastsForExperiment.pl -exp E-MTAB-1066  -conf analysis/differential/reference_assay_group_factor_values.txt

        diffAtlas_generateContrastsForExperiment.pl -exp E-MTAB-1066  -conf analysis/differential/reference_assay_group_factor_values.txt -outdir /homes/kmegy/tmp 


=cut


use strict ;
use Getopt::Long qw(:config no_ignore_case);


## Initialise global $, @ and %
my ($experiment, $conf, $help) ; #arguments
my ($subDirectory, $commandLine) ; #values infered from command line
my $outdir = "./" ;
my %H_config ; #contains info. from the config file 
my %H_hybNameFactorValue ; #store association factor value / hybridization name
my @A_groups ; #store groups 

my $errorCode = 0 ; #report error when testing the replicates, reference etc. 
my $flag ; #to parse config file
my $reference ; #reference factor value to calculate D.E. against


## Get arguments
################
GetOptions( 'help|Help|h|H' => \$help,
	    'exp=s'=> \$experiment,
            'conf=s' => \$conf,
            'out=s'  => \$outdir,
          ) ;

$commandLine = join(' ',@ARGV); 

if (!$experiment || !$conf || !$outdir) { print "[WARNING] Missing experiment (-exp $experiment), configuration files (-conf $conf)\n" ; $help  = 1 ; }

if (!$outdir) { print "[WARNING] No output directory provided (-out $outdir). Default is current directory.\n" } ;

if ($help) { usage($commandLine) ; die ; }


## Directory with SDRF file
my $experimentDirectory = "/net/isilon5/ma/home/arrayexpress/ae2_production/data/EXPERIMENT/" ;

## Experiment sub-directory
if ($experiment =~ /E-(\w+?)-\d+?/) {$subDirectory = $1 ; }
else { die "[ERROR] Experiment -$experiment- not formatted as expected.\n" ; }

## SDRF (input) and XML (output) files
my $sdrf = "$experimentDirectory/$subDirectory/$experiment/$experiment.sdrf.txt" ;
my $outfileXML = "$outdir/$experiment-configuration.xml" ;


## [KM] If I run that program with "$outdir" containing relative path (~/) it works, 
# [KM] ... but it might not work on LSF. 
# [KM] I can easily plan for that (replacing ~/ by $HOME/) - but is it worth doing it?


## Extract information from the config file
###########################################
open (CONF, $conf) || die "Can't open the configuration file $conf!\n" ;
while (my $line=<CONF>) {
        chomp $line ; 

	if ($line !~ /^#/) {
		if ($line =~ /REFERENCE/) { $flag = "REFERENCE" ; }
		elsif ($line =~ /FACTOR_VALUE_KILL/) { $flag = "FACTOR_VALUE_KILL" ; }
		elsif ($line =~ /\[\//) { $flag = "" ; }
        	else { if (($flag ne "") && ($line ne "")) { $H_config{$flag}{$line} = 1 ; } }
	}
}	
close CONF ;


##Print for a test
#foreach my $category (keys %H_Config) {
#	print ">>$category<<\n" ; 
#	foreach my $value (keys %{$H_Config{$category}}) { print ".$value.\n" ; }	
#} print "=====\n" ;
#exit ;



## Collect FactorValues & ENA IDs
## ... to start with, use a basic gawk command
## .... and focuse on a single experiment 
##############################################
#
#Quick command lines
#Later: use Maria's subroutine 
#
my $hybridizationNameFactorValue = `awk -F"\t" '{print \$15"\t"\$NF}' $sdrf` ;
if ($hybridizationNameFactorValue eq "") { die "Couldn't parse SDRF file $sdrf\n" ; } 
#print $hybridizationNameFactorValue ;


#Get hybridization name and factor value from the output of the above command
#Store in a hash of array:
#	$H_Name_FactorValue{FactoreValueA} => [ HybName1, HybName2, HybName3]
#	$H_Name_FactorValue{FactoreValueB} => [ HybName4, HybName5, HybName6, HybName7]
#	etc. 	
my @A_hybridizationNameFactorValue = split ("\n", $hybridizationNameFactorValue) ;	
foreach my $hybridizationNameFactorValue (@A_hybridizationNameFactorValue) {

	my ($hybridizationName, $factorValue) = split ("\t", $hybridizationNameFactorValue) ;
	if ($hybridizationName ne "" && $factorValue ne "" && $hybridizationName !~ "Hybridization") { #discard header or empty lines
		print "Read: -$factorValue- and -$hybridizationName-\n" ;
		push(@{$H_hybNameFactorValue{$factorValue}}, $hybridizationName) ;

	}	
}

## Factor Value parsing
## Parse the factor values: 
#	- >= 3 replicates	  -- delete the Factor Value if not true 
#	- reference		  -- from the list of references generated from the config file
#	- forbidden Factor Values -- from the kill list generated from the config file; delete Factor Value if true

foreach my $FV (keys %H_hybNameFactorValue) {
print "Testing $FV -- @{$H_hybNameFactorValue{$FV}}\n" ;

#Test for forbidden factor value (e.g. 'individual')
if (exists $H_config{"FACTOR_VALUE_KILL"}{$FV}) { delete $H_hybNameFactorValue{$FV} ; print "\tKill List!\n" ; } 		

#Test for reference
if (exists $H_config{"REFERENCE"}{$FV}) { $reference = $FV ; print "\tReference!\n" ; }  

#Test for replicates
my $replicateCount = scalar($H_hybNameFactorValue{$FV}) ;
if ($replicateCount < 3) { delete $H_hybNameFactorValue{$FV} ; print "\tLess than 3 replicates\n" ; }
}

#Anything left afterwards
print "[INFO] Checking Factor Values suitability for differential expression analysis\n" ;

#Reference Factor Value ? 
if (!defined $reference) { print "[INFO][CHECK] Factor Value check FAILED! No reference ($reference)\n" ; $errorCode = 1 ; }

#Any Factor Value left (>= 3 replicates)?
#Need at least 2 of them!
if (keys %H_hybNameFactorValue < 2) { print "[INFO][CHECK] Factor Value check FAILED! Less than 2 values with at least 3 replicates!\n" ; $errorCode = 1 ; }

#If no error reported, then FactorValue test passed!
if ($errorCode == 0) { 

	print "[INFO][CHECK] Factor Value test passed!\n" ; 

	##Make groups 
	# Easier than making them on the fly when generating the XML
	# For easiness, g1 will always be the reference and the rest assigned at random
	my $groupCounter = 1 ;
	$A_groups[$groupCounter] = $reference ; #$A_groups[0] empty, so that array position serves as group ID (g1, g2, g3 etc.) 

	foreach my $FV (keys %H_hybNameFactorValue) {

		if ($FV ne $reference) {
			$groupCounter++ ;
			$A_groups[$groupCounter] = $FV ;
		}		
	}

	##Format in XML
	open (XML, ">$outfileXML") || die ("Can't open output XML file $outfileXML\n") ;
	
	#Beginning XML
	&XMLboundaries("start") ;

	#Assay group section
	&printAssayGroup ;

	#Constrast section
	&printContrast ;

	#End XML
	&XMLboundaries("end") ;

	close XML ;

} else {  #Cannot generate contrast file
	die "[INFO] Contrast file cannot be generated\n" ; 
} 	


## Subroutine
#############
#Print usage for the program
sub usage {
	my ($command_line) = @_;

    	print "Your command line was:\t".
	    "$0 $command_line\n".
	    "Compulsory parameters:\n".
	    "\t-exp: experiment name. E.g. E-MTAB-123\n".
	    "\t-conf: configuration file name\n".
	    "Optional parameters:\n".
	    "\t-outdir: output directory. Default is the current directory.\n" ;
}


#Print tabulations n times, in XML output file 
#'n' being given in argument
sub tabulationXML {
	print XML "\t" x $_[0] ;
}


#Print 'assay_group' section
# [KM] Not sure how to format this to make it easy to read!
# [KM] &tabulationXML(x) same line as the print 'xxx' ???
sub printAssayGroup {
	&tabulationXML(2) ; print XML "<assay_groups>\n" ;
	foreach my $i (1..$#A_groups) { #there is nothing in $A_groups[0]
        	my $factVal = $A_groups[$i] ;

        	&tabulationXML(3) ; print XML "<assay_group id=\"g$i\">\n" ;
        	foreach my $Names (@{$H_hybNameFactorValue{$factVal}}) {
                	&tabulationXML(4) ;
                	print XML "<assay>$Names</assay>\n" ;
        	}
        	&tabulationXML(3) ; print XML "</assay_group>\n" ;

	}
	&tabulationXML(2) ; print XML "</assay_groups>\n" ;
}



#Print 'contrast' section
# [KM] Not sure how to format this to make it easy to read!
# [KM] &tabulationXML(x) same line as the print 'xxx' ???
sub printContrast {

	&tabulationXML(2) ; print XML "<contrasts>\n" ;
	foreach my $i (2..$#A_groups) { #starting at 2 because we have g vs. the rest

        	&tabulationXML(3) ; print XML "<contrast id=\"g1_g".$i."\">\n" ;
        	&tabulationXML(4) ; print XML "<name>'$A_groups[$i]' vs '$A_groups[1]'</name>\n" ;
        	&tabulationXML(4) ; print XML "<reference_assay_group>g1</reference_assay_group>\n" ;
        	&tabulationXML(4) ; print XML "<test_assay_group>g".$i."</test_assay_group>\n" ;
        	&tabulationXML(3) ; print XML "</contrast>\n" ;
	}	
	&tabulationXML(2) ; print XML "</contrasts>\n" ;
}	


#Print the beginning/end of the XML file
#Both are in the same subroutine to make it easier 
#to check that what's been open is being close
sub XMLboundaries {
	my $location = $_[0] ;

	if ($location =~ /start/) {
		print XML "<configuration>\n" ;
		&tabulationXML(1) ; print XML "<analytics>\n" ;
	}	

	if ($location =~ /end/) {
		&tabulationXML(1) ; print XML "</analytics>\n" ;
		print XML "</configuration>\n" ;
	}	
}

