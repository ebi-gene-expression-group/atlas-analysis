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

  none

=head1 EXAMPLES

  diffAtlas_generateContrastsForExperiment.pl -exp experiment_name -conf config_file.txt -out path/to/output/XML/file

  E.g. 
	diffAtlas_generateContrastsForExperiment.pl -exp E-MTAB-1066  -conf analysis/differential/reference_assay_group_factor_values.txt
     or 
        diffAtlas_generateContrastsForExperiment.pl -exp E-MTAB-1066  -conf analysis/differential/reference_assay_group_factor_values.txt -outdir /homes/kmegy/tmp 


=cut


use strict ;
use Getopt::Long qw(:config no_ignore_case);


## Initialise global $, @ and %
my ($experiment, $conf, $help) ; #arguments
my ($subDirectory, $commandLine) ; #values infered from command line
my $outdir = "./" ; #default output directory
my %H_config ; #contains info. from the config file 
my %H_hybNameFactorValue ; #store association factor value / hybridization name
my @A_assayGroups ; #store assay groups 

my $errorCode = 0 ; #report error when testing the replicates, reference etc. 
my $errorMessage ; #error message when testing the replicates, reference etc.
my $flag ; #to parse config file
my $reference ; #reference Factor Value(s) to calculate D.E. against


## Get arguments
################
GetOptions( 'help|Help|h|H' => \$help,
	    'exp=s'=> \$experiment,
            'conf=s' => \$conf,
            'out=s'  => \$outdir,
          ) ;

$commandLine = join(' ',@ARGV); 

if (!$experiment || !$conf) { print "[WARNING] Missing experiment (-exp $experiment) and/or configuration files (-conf $conf)\n" ; $help  = 1 ; }

if (!$outdir) { print "[WARNING] No output directory provided (-out $outdir). Default is current directory.\n" } ;

if ($help) { usage($commandLine) ; die ; }


## Directory with SDRF file
my $experimentDirectory = "/nfs/ma/home/arrayexpress/ae2_production/data/EXPERIMENT/" ;

## Experiment sub-directory
if ($experiment =~ /E-(\w+?)-\d+?/) {$subDirectory = $1 ; }
else { die "[ERROR] Experiment -$experiment- not formatted as expected.\n" ; }

## SDRF (input) and XML (output) files
my $sdrf = "$experimentDirectory/$subDirectory/$experiment/$experiment.sdrf.txt" ;
my $outfileXML = "$outdir/$experiment-configuration.xml.auto" ;


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


##### REMOVE ONCE PROGRAM FINISHED ###
### Print for a test
#foreach my $category (keys %H_Config) {
#	print ">>$category<<\n" ; 
#	foreach my $value (keys %{$H_Config{$category}}) { print ".$value.\n" ; }	
#} print "=====\n" ;
#exit ;
##### END of 'REMOVE ONCE PROGRAM FINISHED'


## Collect FactorValues & ENA IDs
## ... to start with, use Maria's subroutine
## ... later: use MAGETAB module
## ... on a single experiment
##############################################

# Using Maria's subroutine
my ($factorvalueType, $Href_efvs2runAccessions) = &readSDRF($sdrf) ;
my %H_eFactorValues2runIDs = %$Href_efvs2runAccessions ; #dereference the hash

##Print for a test
print "\nProcessing experiment: $experiment..." ;
print "=====> Print for a test\n" ;
foreach my $species (keys %H_eFactorValues2runIDs) {
	print "Species is $species ($factorvalueType)\n" ;
	foreach my $array (keys %{$H_eFactorValues2runIDs{$species}}) {
		print "\tArray is $array\n" ;
		foreach my $factValue (keys %{$H_eFactorValues2runIDs{$species}{$array}}) {
                	print "[REAL] $species - $array - $factValue @{$H_eFactorValues2runIDs{$species}{$array}{$factValue}}\n" ;
		}
	}
}
print "\n\n\n" ;


## Factor Value parsing
## For each organism and each array design (usually: 1 of each only),
## Parse the factor values: 
#	- >= 3 replicates	  -- delete the Factor Value if not true 
#	- reference		  -- from the list of references generated from the config file
#	- forbidden Factor Values -- from the kill list generated from the config file; delete Factor Value if true


#Open the output XML file
open (XML, ">$outfileXML") || die ("Can't open output XML file $outfileXML\n") ;
my $configurationTag = 0 ;
foreach my $species (keys %H_eFactorValues2runIDs) {
        print "Species is $species ($factorvalueType)\n" ;
        foreach my $array (keys %{$H_eFactorValues2runIDs{$species}}) {
                print "\tArray is $array\n" ;

		foreach my $FV (keys %{$H_eFactorValues2runIDs{$species}{$array}}) {
			print "Testing $FV -- @{$H_eFactorValues2runIDs{$species}{$array}{$FV}}\n" ; ##REMOVE ONCE PROGRAM FINISHED

			#Test for forbidden factor value (e.g. 'individual')
			if (exists $H_config{"FACTOR_VALUE_KILL"}{$FV}) { delete $H_eFactorValues2runIDs{$species}{$array}{$FV} ; next ; } 		

			#Test for reference
			if (exists $H_config{"REFERENCE"}{$FV}) { $reference = $FV ; print"\tReference!\n" ; } ##REMOVE print ONCE PROGRAM FINISHED 

			#Test for replicates
			my $replicateCount = scalar @{$H_eFactorValues2runIDs{$species}{$array}{$FV}} ;
			print "Replicate number: $replicateCount\n" ; #FOR TESTING PURPOSE ONLY
			if ($replicateCount < 3) { delete $H_eFactorValues2runIDs{$species}{$array}{$FV} ; print "\tLess than 3 replicates\n" ; } ##REMOVE print ONCE PROGRAM FINISHED
		}	


		#Anything left afterwards
		print "[INFO] Checking Factor Values suitability for differential expression analysis\n" ;

		#Reference Factor Value ? 
		if (!defined $reference) { 
			$errorCode = 1 ;
			$errorMessage .= "No reference. \n" ;
		}

		#Any Factor Value left (>= 3 replicates)?
		#Need at least 2 of them!
		
		if (keys %{$H_eFactorValues2runIDs{$species}{$array}} < 2) {
		#if (keys %H_hybNameFactorValue < 2) {
			$errorCode = 1 ;
			$errorMessage .= "Less than 2 values with at least 3 replicates. \n" ; 
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

			print "=====> Print XML contrast file $outfileXML\n" ;
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
			die "[ERROR] Contrast file cannot be generated for $experiment :: $species :: $array: $errorMessage\n" ; 
		}
	}		 	
}

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
	    "\t-conf: configuration file name\n".
	    "Optional parameters:\n".
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

	$factorvalueType = lc($factorvalueType) ;

	&tabulationXML(2) ; print XML "<contrasts>\n" ;
	foreach my $i (2..$#A_assayGroups) { #starting at 2 because we have g1 (reference) vs. the rest

        	&tabulationXML(3) ; print XML "<contrast id=\"g1_g".$i."\">\n" ;
        	&tabulationXML(4) ; print XML "<name>$factorvalueType:'$A_assayGroups[$i]' vs '$A_assayGroups[1]' on $subArray</name>\n" ;
        	&tabulationXML(4) ; print XML "<reference_assay_group>g1</reference_assay_group>\n" ;
        	&tabulationXML(4) ; print XML "<test_assay_group>g".$i."</test_assay_group>\n" ;
        	&tabulationXML(3) ; print XML "</contrast>\n" ;
	}	
	&tabulationXML(2) ; print XML "</contrasts>\n" ;
}	


#Print the beginning/end of the XML file
#Both are in the same subroutine to make it easier 
#to check that what's been open is being closer
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



#sub readSDRF - modified from the program gxa_summarizeFPKMs.pl by Maria Keays (mkeays@ebi.ac.uk)
# If we are passed the "-a" option, that means we want to calculate the average
# FPKM for each group of replicates (e.g. for each tissue) in the dataset, and
# not just return the FPKM for every run individually.
# To do this we need to find out which run accessions (and hence FPKM tracking files) are
# associated with which factor values. The mapping between run accessions and
# factor values is available in the experiment's SDRF file. Run accessions are
# in a column called "Comment[ENA_RUN]" and factor values in the
# "FactorValue[xxx]" column(s).
#
# Here we read the SDRF file and create a hash mapping run accessions to factor
# values. This must also be done "per organism", in cases where an experiment
# contains data from more than one organism. E.g. we could have data from human
# and mouse liver, but potentially the organism names might not be in a factor
# value column (though they should be?). In this case we would have human ind
# mouse runs assigned "liver" as the factor value and no way to tell them
# apart. So we'll create a key for each organism and after that a key for each
# factor value assigned to each organism. Assigned to each factor value will be
# an array of the run accessions that were found to correspond to them in the
# SDRF.
#
# The (anonymous) hash ends up like e.g.:
#       $efvs2runAccessions->{ "Homo sapiens" }->{ "liver" } = [ "SRR00001", "SRR00002", "SRR00003", ... ]
#                                              ->{ "heart" } = [ "SRR00004", "SRR00005", "SRR00006", ... ]
#
# Arguments:
#       - $logfileHandle : to write to log file.
#       - $sdrfFile : filename of SDRF.
#       - $runAccessions : arrayref of run accessions.
#
# Returns:
#       $efvs2runAccessions : reference to anonymous hash with mappings between factor
#       values (for each organism) and run accessions.
sub readSDRF {
        # Get the arguments:
        my $sdrfFile = $_[0] ; 
       
	#Declare variables
	my ($ENArunInd, $assayNameInd, $arrayNameInd, $arrayDesignInd, $technologyTypeInd) ; #indexes 
	my ($ENArun, $organism, $runAccession, $technologyType, $efvType) ;

        # Simple way to check we have an SDRF.
       if($sdrfFile !~ /\.sdrf\.txt$/) {
                die "$sdrfFile doesn't look like an SDRF file to me. If it is, please append \".sdrf.txt\" to its name.\n";
        }

        # Open SDRF
        open(my $sdrfHandle, "<", $sdrfFile) or die "Can't open $sdrfFile: $!\n";
        
        # Variables to store run name and factor value column indices.
        my ($SDRFrunInd, @efvIndArray, $SDRForgInd);
        
        # Hash reference to remember which run accessions belong to which combination of factor values.
        # $efvs2runAccessions->{ "Gallus gallus" }->{"AFFY-35"}->{ "brain" } = [ "SRR0005", "SRR0006", ... ]    
        my $efvs2runAccessions = {};

        # Remember the run accessions we've seen in the SDRF, because in SDRFs with paired-end
        # data, we might have two rows per run!
        my $seenRuns = [];
        
        # Loop through file.
        while(defined (my $line = <$sdrfHandle>))  {
        	# Remove newlines
        	chomp($line);
        
        	# split on tabs
                my @lineSplit = split("\t", $line);
        
               # Get column indices of run names and factor value(s) from the first line with headers.
               # Characteristics[Organism] should be present in every SDRF header line in some form.
               if($line =~ /characteristics\s*\[\s*organism\s*\]/i) {

               		# Index of Characteristics[Organism] column. Technically it should
               		# be in a FactorValue[] column if it varies but just in case we'll
               		# account for possibility of different organisms here.
                        my @orgIndArray = grep $lineSplit[$_] =~ /characteristics\s*\[\s*organism\s*\]/i, 0..$#lineSplit;
                        $SDRForgInd = $orgIndArray[0];
            
                        # Index of FactorValue[xxxx] column(s)
                        @efvIndArray = grep $lineSplit[$_] =~ /factor\s*value\s*\[/i, 0..$#lineSplit;
                      
			# Need to capture the factor value type (FactorValue[TYPE])
			# ... get the 1st one only to start with
			my $efv = $lineSplit[$efvIndArray[0]] ;
			if ($efv =~ /factor\s*value\s*\[(.+?)\]/i) { $efvType = $1 ; } 

			# Assay names - for old SDRF, use ENA_RUN (RNA-Seq) or Hybridization name (for microarray) 
			#		for newer SDRF, there is an "Assay name" tag
			# Expecting to get assayName from a single of those columns, 
			# ... so using the same $assayName variable		
			# Index of run names from ENA_RUN (old SDRF - RNA-Seq only)             
			my @runIndArray = grep $lineSplit[$_] =~ /ENA_RUN/i, 0..$#lineSplit;
			if (defined $runIndArray[0]) { 
				$ENArunInd = $runIndArray[0];
				$assayNameInd = $runIndArray[0];			
			}

			# Hybridization name (old SDRF - microarray only)
			my @runIndArray = grep $lineSplit[$_] =~ /Hybridization Name/i, 0..$#lineSplit;	
			if (defined $runIndArray[0]) { $assayNameInd = $runIndArray[0]; }
			
			# Assay name (newer SDRF - RNA-Seq and microarray)
			my @runIndArray = grep $lineSplit[$_] =~ /Assay Name/i, 0..$#lineSplit;
                        if (defined $runIndArray[0]) { $assayNameInd = $runIndArray[0]; }

			# Array design
			my @runIndArray = grep $lineSplit[$_] =~ /Array Design/i, 0..$#lineSplit;
			$arrayDesignInd = $runIndArray[0];

			# Technology Type - for old one, infer from ENA_RUN/FASTQ_URI (RNA-Seq) or Hyb Name/Array Design (microarray) 
			# 		    for new ones, there is a "Technology Type" tag
			# Array Names REF (old SDRF - microarray)
			my @runIndArray = grep $lineSplit[$_] =~ /Array Name REF/i, 0..$#lineSplit;
                        $arrayNameInd = $runIndArray[0];

			# Technology Type (newer SDRF - RNA-Seq and microarray)
			my @runIndArray = grep $lineSplit[$_] =~ /Technology Type/i, 0..$#lineSplit;
                        $technologyTypeInd = $runIndArray[0];

			#beware that if no Technology Type tag, get info. from other columns
			#but $technologyType is a string rather than the index (position) of that value
			if (!$technologyTypeInd) {
				if ($arrayNameInd || $assayNameInd) { $technologyType = "array" ; }				
				elsif ( $ENArun ) { $technologyType = "seq" ; }
				#else { die "Cannot identify technology! $technologyTypeInd -$arrayDesignInd;$assayNameInd;$ENArun-\n" ; } #be more informative here? 
			}

		#If line is not header, get the values	
		} else {
			# Check we got indices for organism, run accessions and EFVs
			unless (defined($SDRForgInd) && defined($assayNameInd) && (defined($technologyTypeInd) || defined ($technologyType)) && @efvIndArray) { 
                                die("Do not know SDRF columns for Organism, Assay Name, Technology Type or FactorValues.\n");
			}

                        # Get the run accession and factor values from their indices.
                        my $runAccession = $lineSplit[$assayNameInd];
			$runAccession =~ s/^\"// ; #remove potential initial double quote
			$runAccession =~ s/\"$// ; # remove potential final double quote

			# Get the organism
                        my $organism = $lineSplit[$SDRForgInd];

			# Get the technology type, if not already defined
			if (!defined $technologyType) { $technologyType = $lineSplit[$technologyTypeInd] ; }	

			# If technology is microarray, 
			# Get the array design, otherwise set to 0
			my $arrayDesign = 0 ;
			if ($technologyType =~ /array/) { $arrayDesign = $lineSplit[$arrayDesignInd] ; }
	
                        # Skip to next line if we have already seen this run accession -- SDRFs for
                        # paired-end data have run accessions in twice.
                        if( grep $_ eq $runAccession, @{ $seenRuns } ) { next; }
                        # Otherwise add this run's accession to the seenRuns array so we'll skip it next time.
                        else { push @{ $seenRuns }, $runAccession; }

			print "==> $runAccession <==\n" ; 

                        my @efvArray = @lineSplit[@efvIndArray];
                        my $efvString = join " ", @efvArray;

                        # Add run accession to the right array in %efvs2runAccessions
                        if(exists($efvs2runAccessions->{ $organism }->{ $arrayDesign }->{ $efvString })) {
        			##print "Adding in table $efvs2runAccessions -> $organism -> $efvString , $runAccession\n" ;
	                        push @{ $efvs2runAccessions->{ $organism }->{ $arrayDesign }->{ $efvString } }, $runAccession;
                        }
                        else {
				##print "Adding in table $efvs2runAccessions -> $organism -> $efvString = [ $runAccession ]\n" ;
                                $efvs2runAccessions->{ $organism }->{ $arrayDesign }->{ $efvString } = [ $runAccession ];
                        }
                }
        }
        # Close the file handle
        close($sdrfHandle);


	##Print for testing - to be removed after
	#foreach my $species (keys %$efvs2runAccessions) {
        #	print "[SUB] Species is $species\n" ;
        #	foreach my $array (keys %{$efvs2runAccessions->{$species}}) {
	#		print "[SUB] \tArray is $array\n" ;
        #       	foreach my $organ (keys %{$efvs2runAccessions->{$species}->{$array}}) {
        #        		print "[SUB] $species - $array - $organ %$efvs2runAccessions->{$species}->{$array}->{$organ}\n" ;
	#		}
        #	}
	#}	

       # Return the factor value type (e.g. "genotype") and the reference to a hash of mappings between factor values and run accessions
       return($efvType, $efvs2runAccessions);
}











  
                           

