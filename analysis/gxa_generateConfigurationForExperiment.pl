#!/usr/bin/env perl

# POD documentation - main docs before the code
=pod

=head1 NAME

	Karyn Megy - 18-June-13
	kmegy@ebi.ac.uk
	diffAtlas_generateContrastsForExperiment.pl

=head1 SYNOPSIS

	Generate the XML config file from a IDF/SDRF file

=head1 DESCRIPTION

	Generate the XML config file from a IDF/SDRF file and a generic config file. 

	The generic config file contains the list of reference terms and factor values to exclude from further analysis.
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

	See help for list of options: 
	gxa_generateConfigurationForExperiment.pl -h

=head1 EXAMPLES

	gxa_generateConfigurationForExperiment.pl -exp experiment_name -conf config_file.txt -out path/to/output/XML/dir  -differential|-baseline

	E.g. 
	   gxa_generateConfigurationForExperiment.pl -exp E-MTAB-1066  -conf analysis/differential/reference_assay_group_factor_values.txt -differential -pese PE
     	or 
       gxa_generateConfigurationForExperiment.pl -exp E-MTAB-1066  -conf analysis/differential/reference_assay_group_factor_values.txt -outdir /homes/kmegy/tmp -differential 

=cut


use strict ;
use Getopt::Long qw(:config no_ignore_case) ;
use Magetab4Atlas ;
use XML::Writer ;
use IO::File ;
use EBI::FGPT::Resource::Database; #Atlas DB access


## Initialise global $, @ and %
my ($experiment, $conf, $referenceArg, $killFactorValue, $help, $debug, $differential, $baseline, $noreplicate, $pese) ; #arguments
my ($subDirectory, $commandLine) ; #values infered from command line
my $outdir = "." ; #default output directory
my $experimentType ; #experiment (atlas) type
my %H_config ; #contains info. from the config file 
my @A_assayGroups ; #store assay groups
$A_assayGroups[0] = "none" ; #placeholder, we want assay_groups to start at 1, not 0 


my $flag ; #to parse config file
my $noReferenceError ; # factor: candidate values - to output in log when reporting that no reference could be found     

#to use the XML:writer module
my $XML ;
my $writer ;


## Get arguments
################
GetOptions( 
	'help|Help|h|H' => \$help,
	'exp=s'	 		=> \$experiment,
	'conf=s' 		=> \$conf,
	'differential'	=> \$differential,
	'baseline'		=> \$baseline,
	'noreplicate'	=> \$noreplicate,
	'pese=s'		=> \$pese,
	'ref=s' 		=> \$referenceArg,
	'kill=s' 		=> \$killFactorValue,
	'outdir=s'  	=> \$outdir,
	'debug'			=> \$debug,
) ;

$commandLine = join(' ',@ARGV); 

#Print experiment
##Target error messages when running lots of experiments 
print "[INFO] Generating XML config file for $experiment\n" ;

if (!$experiment || (!$differential && !$baseline) ) { print "[WARNING] Missing experiment (-exp $experiment) or analysis type (-differential or -baseline)\n" ; $help  = 1 ; }
if ($differential && !$conf) { print "[WARNING] Missing configuration files (-conf $conf) for differential analysis\n" ; $help  = 1 ; }
if (!$outdir) { print "[WARNING] No output directory provided (-out $outdir). Default is current directory.\n" } ;
$pese = uc($pese) ;
if ($pese && ($pese ne "PE" && $pese ne "SE")) { print "[WARNING] Value for -pese whould be PE (to restrict to pair end libraries) or SE (to restrict to single end libraries).Value entered: $pese.\n" ; $help = 1 ; }
if ($help) { usage($commandLine) ; die ; }


## Experiment sub-directory
if ($experiment =~ /E-(\w+?)-\d+?/) { $subDirectory = $1 ; }
else { die "[ERROR] Experiment $experiment: name not formatted as expected. Expecting format E-xxx-0123\n" ; }

## Directory with SDRF file
my $experimentDirectory = "/nfs/ma/home/arrayexpress/ae2_production/data/EXPERIMENT/" ;

## IDF & SDRF files (input) and XML file (output)
my $idf = "$experimentDirectory/$subDirectory/$experiment/$experiment.idf.txt" ;
my $outfileXML = "$outdir/$experiment-configuration.xml.auto" ;

## Get list of miRNA
my @A_miRnaList = glob("/ebi/microarray/home/atlas3-production/bioentity_properties/mirbase/*.A-*.tsv") ;
my %H_miRnaList ;

foreach my $miRNA (@A_miRnaList) {
	(my $arrayDesign = $miRNA) =~ s/.*(A-\w{4}-\d+)\.tsv/$1/ ;
	$arrayDesign =~ s/^A-/E-/ ;
	$H_miRnaList{$arrayDesign} = 1 ;
}


## Get array names from Atlas database.
## Extract array ID <-> array name
#########################################################
## Set up database connection
my $dsn = "DBI:Oracle:host=ned.ebi.ac.uk;sid=ATLASPUB;port=1531";
my $username = "atlasprd3";
my $password = "atlas";

# Create connection
if ($debug) { print "[DEBUG] Connecting to Atlas database...\n" } ;
my $atlasDB = EBI::FGPT::Resource::Database->new(
	'dsn' => $dsn,
	'username' => $username,
	'password' => $password,
) or die "Could not connect to Atlas database: $DBI::errstr\n";
if ($debug) { print "[DEBUG] Connected OK.\n"; }

# Get database handle to connect
my $atlasDBH = $atlasDB->get_dbh;

# Create statement handle with query
my $atlasSH = $atlasDBH->prepare("select ACCESSION,NAME from ARRAYDESIGN")
or die "Could not prepare query: ", $atlasDBH->errstr, "\n";

# Execute query
if ($debug) { print "[DEBUG] Querying Atlas database...\n"; }
$atlasSH->execute or die "Could not execute query: ", $atlasSH->errstr, "\n";
if ($debug) { print "[DEBUG] Query successful.\n"; }

# Build hash of results from DB
my %H_arrayIDs2arrayNames ;

# Get each row of results as an arrayref...
while (my $row = $atlasSH->fetchrow_arrayref) {
                                
	# Get the array ID and name
	my ($arrayID, $arrayName) = @{ $row };
                                        
	# And store - complain if already
	if (!exists $H_arrayIDs2arrayNames{$arrayID}) {
		$H_arrayIDs2arrayNames{$arrayID}  = $arrayName ;
	} else {
		die "[ERROR] More than one name for array $arrayID\n" ;
	}
}
# Disconnect from Atlas DB.
$atlasDBH->disconnect;


## For differential analysis only:
## Extract information from the config file
###########################################
if ($differential) {
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
}	

close CONF ;


## Collect experiment_type and ENA IDs, and FactorValues for Differential analysis
# from SDRF file. Use Magetab4atlas module
#################################
# Using readMagetab
my ($expmtType, $factorType, $Href_efvs2runAccessions, $Href_factorValue2factorType, $Href_TechRepsGroup) = &readMagetab($idf) ;

# Dereference hashes
my %H_eFactorValues2runIDs = %$Href_efvs2runAccessions ; 
my %H_factorValue2factorType = %$Href_factorValue2factorType ;
my %H_TechRepsGroup = %$Href_TechRepsGroup ;


# Experiment type
#
# Magetab4atlas module, "get_experiment_type" returns
# Array: "one-colour array" or "two-colour array"
# RNA-Seq: "RNA-seq"
#
# Set to:
#   - microarray_1colour_mrna_differential
#   - microarray_2colour_mrna_differential
#   - microarray_1colour_microrna_differential
#   - rnaseq_mrna_differential
#   - rnaseq_mrna_baseline
my $type ;
if ($expmtType =~ /array/) { $type = "microarray" ;}
if ($expmtType =~ /RNA-seq/) { $type = "rnaseq" ;}

#since color is only for microarray, add the "_", so can be skipped if RNA_seq
my $color ;
if ($expmtType =~ /one-colour array/) { $color = "1colour_" ;}
if ($expmtType =~ /two-colour array/) { $color = "2colour_" ;} 

my $RNA = "mrna";
if (exists $H_miRnaList{$experiment}) { $RNA = "microrna" ;}

my $analysisType ;
if ($differential) { $analysisType = "differential" ;} 
if ($baseline) { $analysisType = "baseline" ;}

#At the moment, we exclude microarray baseline experiments, so die if we have one of those
if ( $type eq "microarray" && $baseline) { die "[ERROR] $experiment - $type experiment cannot be analysed as baseline. Refused for now." ;}

#For each experiment type, make sure I've got the appropritate bits
if ($type =~ /rnaseq/ && $RNA ne '' && $analysisType ne '') {
	$experimentType = "${type}_${color}${RNA}_${analysisType}" ;
	if ($debug) { print "Exepriment (Atlas) type is $experimentType (from \"$expmtType\")\n" ;}
} elsif ($type =~ /microarray/ && $RNA ne '' && $color ne '') {
    $experimentType = "${type}_${color}${RNA}_${analysisType}" ;
	if ($debug) { print "Exepriment (Atlas) type is $experimentType (from \"$expmtType\")\n" ;}
} else { die "[ERROR] $experiment - Cannot get the experiment type: type:$type color:$color RNA:$RNA [from $expmtType]\n" ; } 



############################
## Factor Value parsing
## For each organism and each array design (usually: 1 of each only),
## Parse the factor values: 
#	- >= 3 replicates	      -- delete the Factor Value if not true 
#	- reference		          -- from the list of references generated from the config file, or passed as argument
#	- forbidden Factor Values -- from the kill list generated from the config file, or passed as argument; delete Factor Value if true

my $configurationTag = 0 ;
$noReferenceError = "Candidate reference values for $factorType: ";
if ($debug) { print "[DEBUG] Parsing values collected in Magetab module\n" ; }

#Number of arrays
# - RNA_seq experiment: always one
# - microarray experiment: possibly more than one
my $arrayNumber = scalar keys %H_eFactorValues2runIDs ;
if ($debug) { print "[DEBUG] Number of arrays in that experiment (initially): $arrayNumber\n" ; }

my %H_arrayInAnalyticsElement ; #store arrays that are in an analystics block 
my %H_referenceArray ; #store reference for each couple array/organism

if ($debug) { print "\n[DEBUG] ===== PARSING the DATA =====\n" ; }
foreach my $array (sort keys %H_eFactorValues2runIDs) {
	if ($debug) { print "[DEBUG] Array is $array ($factorType)\n" ; }

	#report error & message when testing the replicates, reference etc. 
	my $errorCode = 0 ;
	my $errorMessage ;

	foreach my $organism (keys %{$H_eFactorValues2runIDs{$array}}) {
		if ($debug) { print "[DEBUG]\tSpecies is $organism\n" ; }

		foreach my $FV (keys %{$H_eFactorValues2runIDs{$array}{$organism}}) {
			if ($debug) { print "[DEBUG] Testing '$FV' -- @{$H_eFactorValues2runIDs{$array}{$organism}{$FV}}\n" ; }
		
			#Differential analysis only
			if ($differential) {
				#Remove forbidden factor value (e.g. 'individual')
				if (exists $H_config{"FACTOR_VALUE_KILL"}{lc($FV)}) { delete $H_eFactorValues2runIDs{$array}{$organism}{$FV} ; next ; } 		
				$noReferenceError .= " '$FV' ";

				#Test for reference - if already one, die loudly
				#(case insensitive: lc only)
				if (exists $H_config{"REFERENCE"}{lc($FV)}) {
					if (!exists $H_referenceArray{$array.$organism}) { $H_referenceArray{$array.$organism} = $FV ;}
					else { 
						$errorCode = 1 ;
						$H_referenceArray{$array.$organism} = "ERROR-MULTIREF" ; #this record a reference error for that array/organism
						$errorMessage .= "More than one reference: $H_referenceArray{$array.$organism} and $FV. " ;
					}
				}
			}

			#Test for statistical replicates
			#Differential: at least 3
			#Baseline: at least 2, unless -noreplicate tag (for specific experiments we really want in Atlas)
			my $replicateCount = scalar @{$H_eFactorValues2runIDs{$array}{$organism}{$FV}} ;
			if ($debug) { print "[DEBUG] Replicate number: $replicateCount\n" ; }
			if ($differential) {
				if ($replicateCount < 3) { 
					delete $H_eFactorValues2runIDs{$array}{$organism}{$FV} ; 
					if ($debug) { print "[DEBUG] Delete H_eFactorValues2runIDs{$array}{$organism}{$FV} due to lack of replicates\n" ; }
				}
			}
           	if ($baseline) {
			   if (($replicateCount < 2) && !$noreplicate) { 
				   delete $H_eFactorValues2runIDs{$array}{$organism}{$FV} ; 
				   if ($debug) { print "[DEBUG] Delete H_eFactorValues2runIDs{$array}{$organism}{$FV} due to lack of replicates\n" ; }	
			   }
		   	}
		}

		#For differential only 
		#Test is reference Factor Value found 
		if ($differential) {
			if (!exists $H_referenceArray{$array.$organism}) {
				$errorCode = 1 ;
				$H_referenceArray{$array.$organism} = "ERROR-NOREF" ; #this record a reference error for that array/organism
				$errorMessage .= "No reference: $noReferenceError. " ;
			}
			
			#Any Factor Value left (>= 3 replicates)?
			#Need at least 2 of them!		
			if (keys %{$H_eFactorValues2runIDs{$array}{$organism}} < 2) {
				$errorCode = 1 ;
				$errorMessage .= "Less than 2 values with at least 3 replicates. " ; 
			}
		}

		#Cannot generate contrast file
		#Warn, dont die
		if($errorCode == 1) { 
			print "\n[WARNING] XML configuration file cannot be generated for $experiment :: $organism :: $array: $errorMessage\n" ;
			if ($debug) { print "[DEBUG] Delete $array ; $organism ...\n" ; }
			delete $H_eFactorValues2runIDs{$array}{$organism} ; 
		}
	}
}
if ($debug) { print "[DEBUG] Finish reading Magetab files\n" ; }


############################
### Now go through the hash again, and print
if ($debug) { print "[DEBUG] Printing XML file\n" ; }

### In the hash, delete any array with no organism left
#this because I need the number of arrays to decided what to print
foreach my $arrCheck (keys %H_eFactorValues2runIDs) {
	if (scalar keys %{$H_eFactorValues2runIDs{$arrCheck}} == 0) {
		delete $H_eFactorValues2runIDs{$arrCheck} ;
	}
}

##Number of arrays
#required later for printing, or not. 
my $arrayNumber = scalar keys %H_eFactorValues2runIDs ;
if ($debug) { print "[DEBUG] Number of arrays in that experiment (after parsing): $arrayNumber\n" ; }
my %H_arrayInAnalyticsElement ; #store arrays that are in an analytics block 

if ($debug) { print "\n[DEBUG] ===== PRINTING the XML =====\n" ; }
foreach my $array (sort keys %H_eFactorValues2runIDs) {
	foreach my $organism (keys %{$H_eFactorValues2runIDs{$array}}) {
		if ($debug) { print "[DEBUG]\tSpecies is $organism\n" ; }

		#If it's the first time, open the XML file
		if ($configurationTag == 0) {
			$configurationTag = 1 ;
			
			#Open file to write to
			$XML = IO::File->new(">$outfileXML");
				
			#Use the XML::Writer module to print
			$writer = XML::Writer->new(OUTPUT => $XML, DATA_MODE => 1, DATA_INDENT => 4);
				
			#Begin configuration XML, add experiment type.
			$writer->startTag("configuration", "experimentType" => $experimentType);
		}

		#If not already done, open the <analytics> element for this array
		if (!exists $H_arrayInAnalyticsElement{$array}) {
			$writer->startTag("analytics");
			$H_arrayInAnalyticsElement{$array} = 1 ;
		}

		### Array_design element - only if experiment has an array
		if ($array ne "0") { $writer->dataElement("array_design" => $array) ; } 

		### Assay_groups element, for each assay
		$writer->startTag("assay_groups");
		   									
		##Make groups
		#Easier than making them on the fly when generating the XML
		#
		#For differential only
		#In a given <assay_groups> element, the lowest index will always be the reference and the rest assigned at random
		my $index = 1 ;
		my $groupCounter = scalar(@A_assayGroups) ; #don't want reset for every organism or array
		my $referenceIndex = 0 ;  #index of the reference Factor Value in @A_assayGroups
			
		if ($differential) {
  			$A_assayGroups[$groupCounter] = $H_referenceArray{$array.$organism} ;
			$referenceIndex = $groupCounter ;
			$index = $groupCounter + 1 ; #because reference is already 1st one
			foreach my $FV (keys %{$H_eFactorValues2runIDs{$array}{$organism}}) {
				if ($FV ne $H_referenceArray{$array.$organism}) {
					$groupCounter++ ;
					$A_assayGroups[$groupCounter] = $FV ;
				}
			}
		}

		#For baseline only
		#No reference, so Group order is at random
		if ($baseline) {
			$index = $groupCounter + 1 ;
			$referenceIndex = $groupCounter ; #this is needed later, for printing
			if ($debug) { print "[DEBUG] Baseline - starting index at $index and groupCounter is $groupCounter\n" ; }
			foreach my $FV (keys %{$H_eFactorValues2runIDs{$array}{$organism}}) {
				$A_assayGroups[$groupCounter] = $FV ;
				$groupCounter++ ;
			}
		}
		print "[INFO] Print XML contrast file $outfileXML\n" ;

		##Assay_group element
		foreach my $i ($referenceIndex..$#A_assayGroups) { 
			my $factVal = $A_assayGroups[$i] ;
			if ($debug) { print "[DEBUG] Organism is $organism and array is $array\n"}
			if ($debug) { print "[DEBUG] [$i] F.V. is $factVal\n"}
			
			##Label
			my $label ;
			foreach my $ft (keys %{$H_factorValue2factorType{$factVal}}) {
				#$label .= "$ft: $H_factorValue2factorType{$factVal}{$ft}; " ;	#label="genotype: stwl bam double mutant"
                $label .= "$H_factorValue2factorType{$factVal}{$ft}; " ;  #label="stwl bam double mutant"
			}
			
			##Remove the trailing
			$label =~ s/; $// ;
			
			##Write element
			$writer->startTag("assay_group", "id" => "g$i", "label" => $label) ;
			foreach my $names (@{$H_eFactorValues2runIDs{$array}{$organism}{$factVal}}) {
				if ($names ne '') {
			
					##Technical replicates, if any
					my $tech_rep ;
					if ($H_TechRepsGroup{$names}) { $writer->dataElement("assay" => $names, "technical_replicateid" => $H_TechRepsGroup{$names}) ; }
					else { $writer->dataElement("assay" => $names) ; }
				}
			}
			$writer->endTag("assay_group") ;
		}
		$writer->endTag("assay_groups") ;

		### Contrast element (differential only)
		##Reference can have any ID but: 
		##  - in the contrast ID the reference ID comes first
		##  - in the contrast name the reference name comes last
		if ($differential) {
		$writer->startTag("contrasts") ;
			
			foreach my $i ($index..$#A_assayGroups) { #starting at 2 because we have g1 (reference) vs. the rest
				$writer->startTag("contrast", "id" => "g${referenceIndex}_g$i") ;

				##if 1 array, or RNA_seq
				if ($arrayNumber == 1) { $writer->dataElement("name" => "'$A_assayGroups[$i]' vs '$A_assayGroups[$referenceIndex]'") ; }

				##if > 1 array
				##user friendly name for the array
				else {
					if ($H_arrayIDs2arrayNames{$array} ne '') {	
						$writer->dataElement("name" => "'$A_assayGroups[$i]' vs '$A_assayGroups[$referenceIndex]' on '$H_arrayIDs2arrayNames{$array}'") ;
					} else { die "[ERROR] $experiment - No user-friendly name for array $array\n" ; } 
				}
				$writer->dataElement("reference_assay_group" => "g$referenceIndex") ;
				$writer->dataElement("test_assay_group" => "g$i") ; 
				$writer->endTag("contrast") ;
			}
			$writer->endTag("contrasts") ;
		}
	}

	##If analytics element open for this array, 
	##Close it.
	if ( ($configurationTag == 1) && (exists $H_arrayInAnalyticsElement{$array}) ) {
		$writer->endTag("analytics");
	}
}

if ($configurationTag == 1) {
	##Close the <configuration> tag, and the file
	$writer->endTag("configuration") ;
	$writer->end ;

	##Close file
	$XML->close ;
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
		"\t-conf: generic configuration file for that program. It should be in the same directory, and possible named reference_assay_group_factor_values.txt\n".
        "\t-baseline/-differential: type of analysis to be performed for Atlas\n".
		"Optional parameters:\n".
		"\t-noreplicate: for baseline analysis, allow Factor Values with less than 2 replicates\n".
		"\t-pese: to restrict to pair end (PE) or single end (SE) libraries. Default: no restriction\n".
		"\t-ref: list of possible reference terms to search for. In double quotes and comma separated if multiple. Take precedence over the config file.\n".
		"\t-kill: list of FactorValue terms to discard (to kill). In double quotes and comma separated if multiple. Take precedence over the config file.\n".
		"\t-outdir: output directory. Default is the current directory.\n".
		"\t-debug: debug mode.\n" ;
}


#Read magetab file (SDRF / IDF files)
sub readMagetab {
    my $efvs2runAccessions = {} ;
	my %H_fv2fvType ;
	my $factorTypeString = "" ;
    my %H_technicalReplicateGroup ;

	if ($debug) { print "\n[DEBUG] ===== READING MAGETAB (readMagetab module) =====\n" ; }
	# Create a Magetab4Atlas object. This reads the MAGETAB documents (via the
	# IDF filename pro vided) to get Atlas-relevant information for each assay in
	# the experiment.
	#
	# Catch the error from the magetab4atlas module
	# ... and only print the 1st line
	my $magetab4atlas ;
	eval { $magetab4atlas = Magetab4Atlas->new( "idf_filename" => $idf );  }; 
	if ($@)  { 
		my @A_moduleErrorMessage = split("\n", $@) ;
		die "[ERROR] $experiment -- Error in the Magetab files. $A_moduleErrorMessage[0]\n" ;
	}

	# Test if we find any assays.
	if(!$magetab4atlas->has_assays) { print "[DEBUG] No assay found!\n" ; } #die "No assay for this experiment!\n"; }
	 
	##Experiment type
	if ($debug) { print "[DEBUG] Experiment type: ".$magetab4atlas->get_experiment_type."\n"; }		
	my $expType = $magetab4atlas->get_experiment_type ; 
	
	my @A_magetabAssay = $magetab4atlas->get_assays ;
	if ($debug) { print "[DEBUG] Assays: $A_magetabAssay[0]\n"; }

	##Get the assay, or die
	if (!@{ $magetab4atlas->get_assays }) { die "[ERROR] $experiment - Cannot extract assay: no name or no factor values\n" }

	foreach my $assay4atlas (@{ $magetab4atlas->get_assays }) {
		if ($debug) { print "[DEBUG] Assays found !\n" ; }

		# Get the organism
		my $organism = $assay4atlas->get_organism() ; 
        if ($debug) { print "[DEBUG]\tOrganism:\n\t\t", $assay4atlas->get_organism, " ($organism) \n" ; }

		# Get the assay name
		# Old SDRF: should be ENA_run (RNA-seq) or hybridisation name (microarray)
		# Newer SDRF: assay name field
		#For RNA-seq, also get the number of fastqURI - should be 1 (SE) or 2 (PE)
		my %H_runAccession ;
		if ($debug) { print "[DEBUG] Assay: ", $assay4atlas->get_name, " (",$assay4atlas->get_name(),")\n" ; }	
		if ($expType eq "RNA-seq") {
			if($assay4atlas->has_fastq_uri_set) {
				foreach my $ENArun (keys %{ $assay4atlas->get_fastq_uri_set }) {
					foreach my $fastqUri (@{ $assay4atlas->get_fastq_uri_set->{ $ENArun } }) {
						$H_runAccession{$ENArun}++ ; #this is to record PE (2) or SE (1)
					}
					if ($debug) { print "[DEBUG]\tENA run: $ENArun ($H_runAccession{$ENArun})\n" ; } 
				}
			}
		} else {
        	$H_runAccession{$assay4atlas->get_name()} = 1 ;
		}

		#Technical replicates, if any
		#They are linked to assay, but storing them associated to the run accession
		# ...as they are going to be printed associated to the run accession
		# (same tech. replicate group for all runs of an assay)
		if ($assay4atlas->has_technical_replicate_group) {
			my $technicalReplicateGroup = $assay4atlas->get_technical_replicate_group ;
			foreach my $run (keys %H_runAccession) {
				$technicalReplicateGroup =~ s/group/t/ ;
				$H_technicalReplicateGroup{$run} = $technicalReplicateGroup ;
				if ($debug) { print "[DEBUG] $run storing technical replicate group $technicalReplicateGroup\n" ; }
			}
		}

		# Get the Factor type(s) and values for this assay
		if ($debug) { print "[DEBUG]\tFactors:\n"; }
		my $H_factors = $assay4atlas->get_factors;
		my $factorValueString = "" ; 
		my $factorTypeString = "" ;
		foreach my $factorType (keys %{ $H_factors }) {

            $factorValueString .= $H_factors->{ $factorType }." " ;
			if ($debug) { print "[DEBUG]\t\t$factorType: ", $H_factors->{ $factorType }, "\n"; }

			#If >= 2 factor type, exclude the ones containing "*Time*"		
			unless (keys %{$H_factors} >= 2 && $factorType =~ /Time/) { $factorTypeString .= "$factorType " ; }
		}

		$factorValueString =~ s/\s+?$// ; #remove trailing space
		$factorTypeString =~ s/\s+?$// ; #remove trailing space
		if ($debug) { print "\t\tfactorValue string: $factorValueString\n"; }
        if ($debug) { print "\t\tfactorType string: $factorTypeString\n"; }

		# For each factorValueString, store the factor types and their value
		# key is factorValueString because it needs to link to another hash which key is also factorValueString
		$H_fv2fvType{$factorValueString} = $H_factors ;
		
		# Get array design for a microarray assay
		my $arrayDesign = 0 ;
		if($expType =~ /array/) {
			if($assay4atlas->has_array_design) {
				$arrayDesign = $assay4atlas->get_array_design() ;
				if ($debug) { print "\tArray design:\n\t\t", $assay4atlas->get_array_design, " ($arrayDesign) \n" ; }
			}
		}

		#Store
		#Go through %H_runAccession to store all runAccession	
		#If PE/SE requirement, parse it here
		foreach my $runAcc (keys %H_runAccession) {
			if ( (!$pese) || ($pese eq "PE" && ($H_runAccession{$runAcc} == 2)) || ($pese eq "SE" && ($H_runAccession{$runAcc} == 1)) ) {
				if(exists($efvs2runAccessions->{ $arrayDesign }->{ $organism }->{ $factorValueString })) {
					push @{ $efvs2runAccessions->{ $arrayDesign }->{ $organism }->{ $factorValueString } }, $runAcc ;
				} else {
					$efvs2runAccessions->{ $arrayDesign }->{ $organism }->{ $factorValueString } = [ $runAcc ] ;
				}
			}	
		}
	}

	#Return values:
	# experiment type, factor value types (e.g. "genotype" - if multiple array: same F.V. types because it's the same SDRF file)
	# and the reference to a hash of mappings between factor values and run accession
	return($expType, $factorTypeString, $efvs2runAccessions, \%H_fv2fvType, \%H_technicalReplicateGroup) ;
}
