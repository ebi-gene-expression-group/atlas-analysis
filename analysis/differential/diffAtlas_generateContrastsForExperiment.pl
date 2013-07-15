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

  diffAtlas_generateContrastsForExperiment.pl -sdrf file.sdrf.txt -conf config_file.txt -out contrast_file.xml

  E.g. 
	diffAtlas_generateContrastsForExperiment.pl -sdrf /net/isilon5/ma/home/arrayexpress/ae2_production/data/EXPERIMENT/MTAB/E-MTAB-1066/E-MTAB-1066.sdrf.txt  -conf analysis/differential/reference_assay_group_factor_values.txt



=cut


use strict ;
use Getopt::Long qw(:config no_ignore_case);


## Initialise global $, @ and %
my ($sdrf, $conf, $xml, $help, $command_line) ; #arguments
my %H_config ; #contains info. from the config file 
my %H_hybNameFactorValue ; #store association factor value / hybridization name
my @A_groups ; #store groups 

my $errorCode = 0 ; #report error when testing the replicates, reference etc. 
my $flag ; #to parse config file
my $reference ; #reference factor value to calculate D.E. against


## Get arguments
################
GetOptions( 'help|Help|h|H' => \$help,
	    'sdrf=s' => \$sdrf,
            'xml=s'  => \$xml,	
            'conf=s' => \$conf
          ) ;

$command_line = join(' ',@ARGV); 

if (!$sdrf || !$conf) { print "[WARNING] Missing sdrf (-sdrf $sdrf) or configuration files (-conf $conf)!\n" ; $help  = 1 ; }

#Suppress if automatic output file name generation, or merge with the above warning message
if (!$xml) { print "[WARNING] Missing xml file name (-xml $xml)\n" ; $help = 1 ; }

if ($help) { usage($command_line) ; die ; }

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
	open (XML, ">$xml") || die ("Can't open output XML file $xml\n") ;
	
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
	"\t-sdrf: SDRF file name\n".
	"\t-conf: configuration file name\n".
	"\t-xml: xml output file name\n" ;
}


#Print tabulation n times 
#'n' being given in argument
sub tabulation {
	print XML "\t" x $_[0] ;
}


#Print 'assay_group' section
# [KM] Not sure how to format this to make it easy to read!
# [KM] &tabulation(x) same line as the print 'xxx' ???
sub printAssayGroup {
	&tabulation(2) ; print XML "<assay_groups>\n" ;
	foreach my $i (1..$#A_groups) { #there is nothing in $A_groups[0]
        	my $factVal = $A_groups[$i] ;

        	&tabulation(3) ; print XML "<assay_group id=\"g$i\">\n" ;
        	foreach my $Names (@{$H_hybNameFactorValue{$factVal}}) {
                	&tabulation(4) ;
                	print XML "<assay>$Names</assay>\n" ;
        	}
        	&tabulation(3) ; print XML "</assay_group>\n" ;

	}
	&tabulation(2) ; print XML "</assay_groups>\n" ;
}



#Print 'contrast' section
# [KM] Not sure how to format this to make it easy to read!
# [KM] &tabulation(x) same line as the print 'xxx' ???
sub printContrast {

	&tabulation(2) ; print XML "<contrasts>\n" ;
	foreach my $i (2..$#A_groups) { #starting at 2 because we have g vs. the rest

        	&tabulation(3) ; print XML "<contrast id=\"g1_g".$i."\">\n" ;
        	&tabulation(4) ; print XML "<name>'$A_groups[$i]' vs '$A_groups[1]'</name>\n" ;
        	&tabulation(4) ; print XML "<reference_assay_group>g1</reference_assay_group>\n" ;
        	&tabulation(4) ; print XML "<test_assay_group>g".$i."</test_assay_group>\n" ;
        	&tabulation(3) ; print XML "</contrast>\n" ;
	}	
	&tabulation(2) ; print XML "</contrasts>\n" ;
}	


#Print the beginning/end of the XML file
#Both are in the same subroutine to make it easier 
#to check that what's been open is being close
sub XMLboundaries {
	my $location = $_[0] ;

	if ($location =~ /start/) {
		print XML "<configuration>\n" ;
		&tabulation(1) ; print XML "<analytics>\n" ;
	}	

	if ($location =~ /end/) {
		&tabulation(1) ; print XML "</analytics>\n" ;
		print XML "</configuration>\n" ;
	}	
}

