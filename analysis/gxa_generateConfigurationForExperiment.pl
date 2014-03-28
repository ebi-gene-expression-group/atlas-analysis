#!/usr/bin/env perl

# POD documentation - main docs before the code
=pod

=head1 NAME

	Karyn Megy - 18-June-13
	kmegy@ebi.ac.uk
	gxa_generateConfigurationForExperiment.pl

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

	gxa_generateConfigurationForExperiment.pl -exp expAcc -conf config_file.txt -outdir path/to/output/XML/dir  -differential|-baseline

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
use EBI::FGPT::Resource::Database::GXA; #Atlas DB access
use Cwd;


## Initialise global $, @ and %
my ($expAcc, $conf, $referenceArg, $killFactorValue, $differential, $baseline, $noreplicate, $pese, $help, $debug) ; #arguments
my ($subDirectory, $commandLine) ; #values infered from command line
my $outdir = $expAcc ; #output directory - default is experiment
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
	'exp=s'	 		=> \$expAcc,
	'conf=s' 		=> \$conf,
	'differential'	=> \$differential,
	'baseline'		=> \$baseline,
	'noreplicate'	=> \$noreplicate,
	'pese=s'		=> \$pese,
	'ref=s' 		=> \$referenceArg,
	'kill=s' 		=> \$killFactorValue,
    'out=s'			=> \$outdir,
	'outdir=s'  	=> \$outdir,
	'debug'			=> \$debug,
) ;

$commandLine = join(' ',@ARGV); 

if (!$expAcc) {
  die "[ERROR] Missing experiment (-exp $expAcc)";
}

my $dir = getcwd;
my $logfile = "$dir/atlas_configuration_generation_$expAcc.idf.txt.log" ;
open(my $logFileHandle, ">", $logfile) or die "[ERROR] Can't create log file: $logfile ($!)"; 

if (!$differential && !$baseline)  { 
  &log($logFileHandle, "[ERROR] analysis type (-differential or -baseline)") ; $help  = 1 ; }
if ($differential && !$conf) { 
  &log($logFileHandle, "[ERROR] Missing configuration files (-conf $conf) for differential analysis") ; $help  = 1 ; }
$pese = uc($pese) ;
if ($pese && ($pese ne "PE" && $pese ne "SE")) { 
  &log($logFileHandle, "[ERROR] Value for -pese whould be PE (to restrict to pair end libraries) or SE (to restrict to single end libraries).Value entered: $pese.") ; $help = 1 ; }
if ($help) { usage($commandLine) ; die ; }

&log($logFileHandle, "[INFO] Generating XML config file for $expAcc"); 

## Experiment sub-directory
if ($expAcc =~ /E-(\w+?)-\d+?/) { $subDirectory = $1 ; }
else { &log($logFileHandle, "[ERROR] Experiment $expAcc: name not formatted as expected. Expecting format E-xxx-0123"); die ; }

## Directory with SDRF file
my $experimentDirectory = "/nfs/ma/home/arrayexpress/ae2_production/data/EXPERIMENT/" ;

## IDF & SDRF files (input) and XML file (output)
my $idf = "$experimentDirectory/$subDirectory/$expAcc/$expAcc.idf.txt" ;
my $outfileXML = "$outdir/$expAcc-configuration.xml.auto" ;

## Get list of miRNA
my @A_miRnaList = glob("/ebi/microarray/home/atlas3-production/bioentity_properties/mirbase/*.A-*.tsv") ;
my %H_miRnaList ;

foreach my $miRNA (@A_miRnaList) {
	(my $arrayDesign = $miRNA) =~ s/.*(A-\w{4}-\d+)\.tsv/$1/ ;
	$arrayDesign =~ s/^A-/E-/ ;
	$H_miRnaList{$arrayDesign} = 1 ;
}

## Output directory is the experiment one, by default
# Check if exists, warn, and die, if it doesn't
unless (-d $outdir) {
	&log($logFileHandle, "[ERROR] Output directory ($outdir) doesn't exist. Create it (mkdir $outdir) or give directory name in argument (-out/-outdir)."); die ;
}

## Get array names from Atlas database.
## Extract array ID <-> array name
#########################################################
## Set up database connection
<<<<<<< HEAD
=======
my $dsn = "DBI:Oracle:host=ora-vm-025.ebi.ac.uk;sid=ATLASREL;port=1531";
my $username = "atlasprd3";
my $password = "atlas";
>>>>>>> e885464df4ef4b3ee0b0f765db6048f6bbef44cb

# Create connection
if ($debug) { &log($logFileHandle, "[DEBUG] Connecting to Atlas database...") ; } ;
my $atlasDB = EBI::FGPT::Resource::Database::GXA->new();
if (!$atlasDB) { &log($logFileHandle, "[ERROR] Could not connect to Atlas database: $DBI::errstr"); die ; } ;
if ($debug) { &log($logFileHandle, "[DEBUG] Connected OK.") ; }

# Get database handle to connect
my $atlasDBH = $atlasDB->get_dbh;

# Create statement handle with query
my $atlasSH = $atlasDBH->prepare("select ACCESSION,NAME from ARRAYDESIGN");
if (!$atlasSH) { &log($logFileHandle, "[ERROR] Could not prepare query: ".$atlasDBH->errstr); die ; } ;

# Execute query
if ($debug) { &log($logFileHandle, "[DEBUG] Querying Atlas database...") ; }
if (!$atlasSH->execute) {  &log($logFileHandle, "Could not execute query: ".$atlasSH->errstr); die ; }  ;
if ($debug) { &log($logFileHandle, "[DEBUG] Query successful.") ; }

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
		&log($logFileHandle, "[ERROR] More than one name for array $arrayID"); die ;
	}
}
# Disconnect from Atlas DB.
$atlasDBH->disconnect;


## For differential analysis only:
## Extract information from the config file
###########################################
if ($differential) {
        if (!open (CONF, $conf)) { &log($logFileHandle, "[ERROR] Can't open the configuration file $conf!"); die ; } ;
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
if (exists $H_miRnaList{$expAcc}) { $RNA = "microrna" ;}

my $analysisType ;
if ($differential) { $analysisType = "differential" ;} 
if ($baseline) { $analysisType = "baseline" ;}

#At the moment, we exclude microarray baseline experiments, so die, if we have one of those
if ( $type eq "microarray" && $baseline) { 
  &log($logFileHandle, "[ERROR] $expAcc - $type experiment cannot be analysed as baseline. Refused for now."); die ; }

#For each experiment type, make sure I've got the appropritate bits
if ($type =~ /rnaseq/ && $RNA ne '' && $analysisType ne '') {
	$experimentType = "${type}_${color}${RNA}_${analysisType}" ;
	if ($debug) { &log($logFileHandle, "Exepriment (Atlas) type is $experimentType (from \"$expmtType\")") ;}
} elsif ($type =~ /microarray/ && $RNA ne '' && $color ne '') {
    $experimentType = "${type}_${color}${RNA}_${analysisType}" ;
	if ($debug) { &log($logFileHandle, "Exepriment (Atlas) type is $experimentType (from \"$expmtType\")") ;}
} else { &log($logFileHandle, "[ERROR] $expAcc - Cannot get the experiment type: type:$type color:$color RNA:$RNA [from $expmtType]"); die ; } 



############################
## Factor Value parsing
## For each organism and each array design (usually: 1 of each only),
## Parse the factor values: 
#	- >= 3 replicates	      -- delete the Factor Value if not true 
#	- reference		          -- from the list of references generated from the config file, or passed as argument
#	- forbidden Factor Values -- from the kill list generated from the config file, or passed as argument; delete Factor Value if true

my $configurationTag = 0 ;
$noReferenceError = "Candidate reference values for $factorType: ";
if ($debug) { &log($logFileHandle, "[DEBUG] Parsing values collected in Magetab module") ; }

#Number of arrays
# - RNA_seq experiment: always one
# - microarray experiment: possibly more than one
my $arrayNumber = scalar keys %H_eFactorValues2runIDs ;
if ($debug) { &log($logFileHandle, "[DEBUG] Number of arrays in that experiment (initially): $arrayNumber") ; }

my %H_arrayInAnalyticsElement ; #store arrays that are in an analystics block 
my %H_referenceArray ; #store reference for each couple array/organism

if ($debug) { &log($logFileHandle, "\n[DEBUG] ===== PARSING the DATA =====") ; }
foreach my $array (sort keys %H_eFactorValues2runIDs) {
	if ($debug) { &log($logFileHandle, "[DEBUG] Array is $array ($factorType)") ; }

	#report error when testing the replicates, reference etc. 
	my $warningMessage = "" ;

	foreach my $organism (keys %{$H_eFactorValues2runIDs{$array}}) {
		if ($debug) { &log($logFileHandle, "[DEBUG]\tSpecies is $organism") ; }

		foreach my $FV (keys %{$H_eFactorValues2runIDs{$array}{$organism}}) {
			if ($debug) { &log($logFileHandle, "[DEBUG] Testing '$FV' -- @{$H_eFactorValues2runIDs{$array}{$organism}{$FV}}") ; }
		
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
						$H_referenceArray{$array.$organism} = "ERROR-MULTIREF" ; #this record a reference error for that array/organism
						$warningMessage .= "More than one reference: $H_referenceArray{$array.$organism} and $FV. " ;
					}
				}
			}

			#Test for statistical replicates
			#Differential: at least 3
			#Baseline: at least 2, unless -noreplicate tag (for specific experiments we really want in Atlas)
			my $replicateCount = scalar @{$H_eFactorValues2runIDs{$array}{$organism}{$FV}} ;
			if ($debug) { &log($logFileHandle, "[DEBUG] Replicate number: $replicateCount") ; }
			if ($differential) {
				if ($replicateCount < 3) { 
					delete $H_eFactorValues2runIDs{$array}{$organism}{$FV} ; 
					if ($debug) { &log($logFileHandle, "[DEBUG] Delete H_eFactorValues2runIDs{$array}{$organism}{$FV} due to lack of replicates") ; }
				}
			}
           	if ($baseline) {
			   if (($replicateCount < 2) && !$noreplicate) { 
				   delete $H_eFactorValues2runIDs{$array}{$organism}{$FV} ; 
				   if ($debug) { &log($logFileHandle, "[DEBUG] Delete H_eFactorValues2runIDs{$array}{$organism}{$FV} due to lack of replicates") ; }	
			   }
		   	}
		}

		#For differential only 
		#Test is reference Factor Value found 
		if ($differential) {
			if (!exists $H_referenceArray{$array.$organism}) {
				$H_referenceArray{$array.$organism} = "ERROR-NOREF" ; #this record a reference error for that array/organism
				$warningMessage .= "No reference: $noReferenceError. Add it manually with -ref \"new reference word\"" ;
			}
			
			#Any Factor Value left (>= 3 replicates)?
			#Need at least 2 of them!		
			if (keys %{$H_eFactorValues2runIDs{$array}{$organism}} < 2) {
				$warningMessage .= "Less than 2 values with at least 3 replicates. " ; 
			}
		}

		#Cannot generate contrast file
		#Warn, dont die
		if($warningMessage ne "") { 
			&log($logFileHandle, "\n[WARNING] XML configuration file cannot be generated for $expAcc :: $organism :: $array: $warningMessage") ;
			if ($debug) { &log($logFileHandle, "[DEBUG] Delete $array ; $organism ...") ; }
			delete $H_eFactorValues2runIDs{$array}{$organism} ; 
		}
	}
}
if ($debug) { &log($logFileHandle, "[DEBUG] Finish reading Magetab files") ; }


############################
### Now go through the hash again, and print
if ($debug) { &log($logFileHandle, "[DEBUG] Printing XML file") ; }

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
if ($debug) { &log($logFileHandle, "[DEBUG] Number of arrays in that experiment (after parsing): $arrayNumber") ; }
my %H_arrayInAnalyticsElement ; #store arrays that are in an analytics block 

if ($debug) { &log($logFileHandle, "\n[DEBUG] ===== PRINTING the XML =====") ; }
foreach my $array (sort keys %H_eFactorValues2runIDs) {
	foreach my $organism (keys %{$H_eFactorValues2runIDs{$array}}) {
		if ($debug) { &log($logFileHandle, "[DEBUG]\tSpecies is $organism") ; }

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
		#This will create <assay_groups> for each organisms, for exemple
		#... but this is incorrect!!
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
			if ($debug) { &log($logFileHandle, "[DEBUG] Baseline - starting index at $index and groupCounter is $groupCounter") ; }
			foreach my $FV (keys %{$H_eFactorValues2runIDs{$array}{$organism}}) {
				$A_assayGroups[$groupCounter] = $FV ;
				$groupCounter++ ;
			}
		}
		&log($logFileHandle, "[INFO] Print XML contrast file $outfileXML") ;

		##Assay_group element
		foreach my $i ($referenceIndex..$#A_assayGroups) { 
			my $factVal = $A_assayGroups[$i] ;
			if ($debug) { &log($logFileHandle, "[DEBUG] Organism is $organism and array is $array"); }
			if ($debug) { &log($logFileHandle, "[DEBUG] [$i] F.V. is $factVal"); }
			
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
	
        #This will create <assay_groups> for each organisms, for exemple
		#.. incorrect, but keeping this in was it break somethign else
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
					} else { &log($logFileHandle, "[ERROR] $expAcc - No user-friendly name for array $array"); die ; } 
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

close($logFileHandle);


## Subroutine
#############
#Print usage for the program
sub usage {
	my ($command_line) = @_;
	
	&log($logFileHandle, "Your command line was:\t".
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
		"\t-outdir: output directory. Default is the experiment directory.\n".
        "\t-out: same as -outdir.\n".
		"\t-debug: debug mode.") ;
}


#Read magetab file (SDRF / IDF files)
sub readMagetab {
    my $efvs2runAccessions = {} ;
	my %H_fv2fvType ;
	my $factorTypeString = "" ;
    my %H_technicalReplicateGroup ;

	if ($debug) { &log($logFileHandle, "\n[DEBUG] ===== READING MAGETAB (readMagetab module) =====") ; }
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
		&log($logFileHandle, "[ERROR] $expAcc -- Error in the Magetab files. $A_moduleErrorMessage[0]"); die ;
	}

	# Test if we find any assays.
	if(!$magetab4atlas->has_assays) { &log($logFileHandle, "[DEBUG] No assay found!") ; }
	 
	##Experiment type
	if ($debug) { &log($logFileHandle, "[DEBUG] Experiment type: ".$magetab4atlas->get_experiment_type.""); }		
	my $expType = $magetab4atlas->get_experiment_type ; 
	
	my @A_magetabAssay = $magetab4atlas->get_assays ;
	if ($debug) { &log($logFileHandle, "[DEBUG] Assays: $A_magetabAssay[0]"); }

	##Get the assay, or die
	if (!@{ $magetab4atlas->get_assays }) { &log($logFileHandle, "[ERROR] $expAcc - Cannot extract assay: no name or no factor values"); die ; }

	foreach my $assay4atlas (@{ $magetab4atlas->get_assays }) {
		if ($debug) { &log($logFileHandle, "[DEBUG] Assays found !") ; }

		# Get the organism
		my $organism = $assay4atlas->get_organism() ; 
        if ($debug) { &log($logFileHandle, "[DEBUG]\tOrganism:\n\t\t", $assay4atlas->get_organism, " ($organism) ") ; }

		# Get the assay name
		# Old SDRF: should be ENA_run (RNA-seq) or hybridisation name (microarray)
		# Newer SDRF: assay name field
		#For RNA-seq, also get the number of fastqURI - should be 1 (SE) or 2 (PE)
		my %H_runAccession ;
		if ($debug) { &log($logFileHandle, "[DEBUG] Assay: ", $assay4atlas->get_name, " (",$assay4atlas->get_name(),")") ; }	
		if ($expType eq "RNA-seq") {
			if($assay4atlas->has_fastq_uri_set) {
				foreach my $ENArun (keys %{ $assay4atlas->get_fastq_uri_set }) {
					foreach my $fastqUri (@{ $assay4atlas->get_fastq_uri_set->{ $ENArun } }) {
						$H_runAccession{$ENArun}++ ; #this is to record PE (2) or SE (1)
					}
					if ($debug) { &log($logFileHandle, "[DEBUG]\tENA run: $ENArun ($H_runAccession{$ENArun})") ; } 
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
				$technicalReplicateGroup =~ s/\s//g ;
				$H_technicalReplicateGroup{$run} = $technicalReplicateGroup ;
				if ($debug) { &log($logFileHandle, "[DEBUG] $run storing technical replicate group $technicalReplicateGroup") ; }
			}
		}

		# Get the Factor type(s) and values for this assay
		if ($debug) { &log($logFileHandle, "[DEBUG]\tFactors:"); }
		my $H_factors = $assay4atlas->get_factors;
		my $factorValueString = "" ; 
		my $factorTypeString = "" ;
		foreach my $factorType (keys %{ $H_factors }) {

            $factorValueString .= $H_factors->{ $factorType }." " ;
			if ($debug) { &log($logFileHandle, "[DEBUG]\t\t$factorType: ", $H_factors->{ $factorType }, ""); }

			#If >= 2 factor type, exclude the ones containing "*Time*"		
			unless (keys %{$H_factors} >= 2 && $factorType =~ /Time/) { $factorTypeString .= "$factorType " ; }
		}

		$factorValueString =~ s/\s+?$// ; #remove trailing space
		$factorTypeString =~ s/\s+?$// ; #remove trailing space
		if ($debug) { &log($logFileHandle, "\t\tfactorValue string: $factorValueString"); }
        if ($debug) { &log($logFileHandle, "\t\tfactorType string: $factorTypeString"); }

		# For each factorValueString, store the factor types and their value
		# key is factorValueString because it needs to link to another hash which key is also factorValueString
		$H_fv2fvType{$factorValueString} = $H_factors ;
		
		# Get array design for a microarray assay
		my $arrayDesign = 0 ;
		if($expType =~ /array/) {
			if($assay4atlas->has_array_design) {
				$arrayDesign = $assay4atlas->get_array_design() ;
				if ($debug) { &log($logFileHandle, "\tArray design:\n\t\t", $assay4atlas->get_array_design, " ($arrayDesign) ") ; }
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

# Dual log to STDOUT and $logFileHandle. The latter will be reported by Conan in the interface; the former will 
# be reported in the LSF output file (and via an email from Conan) in the case of failure of this script.
sub log() {
  my $logFileHandle = shift;
  my $msg = shift; 
  print STDOUT "$msg\n";
  printf($logFileHandle "$msg\n");
}

