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

------------------
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

------------------
	There are many hash & arrays in this program;
	The most important ones are: 

	I. Making assay_groups
	 -> main program
	 %H_referenceArray{$array;$organism;$time}{$FVlist} = 1      ... only if F.V. is a reference
	 %H_assayGroups1Difference{$array;$organism;$time}{$FVstrinag1}{$FVstring2} = 1 ... F.V. list with only 1 F.V. difference 
	 %H_inAssayGroup{$array;$organism;$time}{$FVlist} = #        ... F.V. and their assay_group ID (g1, g2 etc.) 

	II. Run_accession/assay_name & their associated factor value 
	 -> from readMagetab
	 %H_eFactorValues2runIDs{$array_design}{$organism}{$time}{$FVstring} = [$runAccession/assay_name]
		(in subroutine: %efvs2runAccess)

	 %H_assayFactorValues{$array_design;$organism;$time}{$assay_name} = $FVstring
		(in subroutine: %H_assayFV)

=head1 OPTIONS

	See help for list of options: 
	gxa_generateConfigurationForExperiment.pl -h

=head1 EXAMPLES

	gxa_generateConfigurationForExperiment.pl -exp expAcc -conf config_file.txt -differential|-baseline

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
use Data::Dumper ;
use Cwd;


## Initialise global $, @ and %
my ($expAcc, $conf, $referenceArg, $killFactorValue, $differential, $baseline, $noreplicate, $pese, $help, $debug) ; #arguments
my ($subDirectory, $commandLine) ; #values infered from command line
my $outdir = $expAcc ; #output directory - default is experiment
my $experimentType ; #experiment (atlas) type
my %H_config ; #contains info. from the config file 
my %H_assayGroups1Difference ; #built the assay_groupS
my %H_assayGroups1DifferenceTimeSeries ; #built the assay_groupS for time series
my %H_inAssayGroup ; #record F.V. and their assay_group
my %H_assayGroupDone ; #store assay_group that have been processed (ALL vs. ALL)
my $differentialTag ; #reference or non-reference
my $flag ; #to parse config file
my $noReferenceError ; # factor: candidate values - to output in log when reporting that no reference could be found     

## Time points with < 3 F.V. sets - specific handling
my %H_timePointsNotEnoughFvSets ; #time points with < 3 F.V. sets 
my %H_timePointsNotEnoughFvSetsFV2runIDs ; #same as %H_eFV2runIDs but limited to time points with < 3 F.V. sets 
my @A_timePointsNotEnoughFvSets ; #order time points 

## For using the XML:writer module
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
	'outdir=s'  	=> \$outdir,
    'out=s'			=> \$outdir,
	'debug'			=> \$debug,
) ;

$commandLine = join(' ',@ARGV); 


## Open the log file
my $dir = getcwd;
my $logfile = "$dir/atlas_configuration_generation_$expAcc.idf.txt.log" ;
open(my $logFileHandle, ">", $logfile) or die "[ERROR] Can't create log file: $logfile ($!)";


## Check arguments, and print any errors in the log file 
if (!$outdir) { 
	&log($logFileHandle, "[WARNING] No output directory provided (-out $outdir). Default is current directory.\n") ; }

if (!$expAcc) {
	&log($logFileHandle, "[ERROR] Missing experiment (-exp $expAcc)") ; $help  = 1 ;  }

if (!$differential && !$baseline)  {
	&log($logFileHandle, "[ERROR] Missing analysis type (-differential or -baseline)") ; $help  = 1 ; }

if ($differential && !$conf) {
	&log($logFileHandle, "[ERROR] Missing configuration files (-conf $conf) for differential analysis") ; $help  = 1 ; }

$pese = uc($pese) ;
if ($pese && ($pese ne "PE" && $pese ne "SE")) {
	&log($logFileHandle, "[ERROR] Value for -pese whould be PE (to restrict to pair end libraries) or SE (to restrict to single end libraries).Value entered: $pese.") ; $help = 1 ; }

if ($help) { usage($commandLine) ; die ; }

&log($logFileHandle, "[INFO] Generating XML config file for $expAcc") ;


## Experiment sub-directory
if ($expAcc =~ /E-(\w+?)-\d+?/) { $subDirectory = $1 ; }
else { &log($logFileHandle, "[ERROR] Experiment $expAcc: name not formatted as expected. Expecting format E-xxx-0123") ; die ; }

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
    &log($logFileHandle, "[ERROR] Output directory ($outdir) doesn't exist. Create it (mkdir $outdir) or give directory name in argument (-out/-outdir).") ; die ;
}


## Get array names from Atlas database.
## Extract array ID <-> array name
#########################################################
## Set up database connection
my $dsn = "DBI:Oracle:host=ora-vm-025.ebi.ac.uk;sid=ATLASREL;port=1531";
my $username = "atlasprd3";
my $password = "atlas";

## Create connection
if ($debug) { &log($logFileHandle, "[DEBUG] Connecting to Atlas database...\n") ; } 
my $atlasDB = EBI::FGPT::Resource::Database->new(
	'dsn' => $dsn,
	'username' => $username,
	'password' => $password,
) ;
if (!$atlasDB) { &log($logFileHandle, "[ERROR] Could not connect to Atlas database: $DBI::errstr") ; die ; }
if ($debug) { &log($logFileHandle, "[DEBUG] Connected OK.") ; }

## Get database handle to connect
my $atlasDBH = $atlasDB->get_dbh ;

## Create statement handle with query
my $atlasSH = $atlasDBH->prepare("select ACCESSION,NAME from ARRAYDESIGN") ;
if (!$atlasSH) { &log($logFileHandle, "[ERROR] Could not prepare query: ".$atlasDBH->errstr) ; die ; }

## Execute query
if ($debug) { &log($logFileHandle, "[DEBUG] Querying Atlas database...") ; }
$atlasSH->execute or die "Could not execute query: ", $atlasSH->errstr, "\n";
if ($debug) { &log($logFileHandle, "[DEBUG] Query successful.") ; }

## Build hash of results from DB
my %H_arrayIDs2arrayNames ;

## Get each row of results as an arrayref...
while (my $row = $atlasSH->fetchrow_arrayref) {
                                
	## Get the array ID and name
	my ($arrayID, $arrayName) = @{ $row };
                                        
	## And store - complain if already
	if (!exists $H_arrayIDs2arrayNames{$arrayID}) {
		$H_arrayIDs2arrayNames{$arrayID}  = $arrayName ;
	} else {
        &log($logFileHandle, "[ERROR] More than one name for array $arrayID") ; die ;
	}
}
## Disconnect from Atlas DB.
$atlasDBH->disconnect;


## For differential analysis only:
## Extract information from the config file
###########################################
if ($differential) {
	if (!open (CONF, $conf)) { &log($logFileHandle, "[ERROR] Can't open the configuration file $conf!") ; die ; } ;
	while (my $line=<CONF>) {
		chomp $line ; 

		if ($line !~ /^#/) {

			## If reference value(s) or kill terms given in argument,
			# they take precedence over the ones from the config file.
			# Should be in double quotes and comma separated
			if ($referenceArg) {
				my @A_referenceArg = split(",", $referenceArg) ; 
				for my $refArg (@A_referenceArg) {
					## For each value, trim any start/end spaces (but not middle ones!)
					# ...and store
					$refArg =~ s/^\s+//g ; $refArg =~ s/\s+$//g ;
					$H_config{"REFERENCE"}{lc($refArg)} = 1 ; 
				}	
			}	
			if ($killFactorValue) { 	
				my @A_killFV = split(",", $killFactorValue) ;
				for my $factorValue (@A_killFV) {
					## For each value, trim any start/end spaces (but not middle ones!)
					# ...and store
					$factorValue =~ s/^\s+//g ; $factorValue =~ s/\s+$//g ;
					$H_config{"FACTOR_VALUE_KILL"}{lc($factorValue)} = 1 ; 
				}
			}

			## If no reference value(s) or kill terms  given in argument,
			# take them from the config file.
			if ($line =~ /REFERENCE/ && !$referenceArg) { $flag = "REFERENCE" ; }
			elsif ($line =~ /FACTOR_VALUE_KILL/ && !$killFactorValue) { $flag = "FACTOR_VALUE_KILL" ; }
			elsif ($line =~ /\[\//) { $flag = "" ; }
			else { if (($flag ne "") && ($line ne "")) { $H_config{$flag}{lc($line)} = 1 ; } } #Use lc for case insensitive comparison
		}
	}	
}	
close CONF ;


## Collect various information from the SDRF file
################################
## Using readMagetab
my ($expmtType, $factorType, $Href_efvs2runAccessions, $Href_factorValue2factorType, $Href_TechRepsGroup, $Href_assayFactorValues) = &readMagetab($idf) ;


## Dereference hashes
my %H_eFactorValues2runIDs = %$Href_efvs2runAccessions ; 
my %H_factorValue2factorType = %$Href_factorValue2factorType ;
my %H_TechRepsGroup = %$Href_TechRepsGroup ;
my %H_assayFactorValues = %$Href_assayFactorValues ;


## Experiment type
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

## Since color is only for microarray, add the "_", so can be skipped if RNA_seq
my $color ;
if ($expmtType =~ /one-colour array/) { $color = "1colour_" ;}
if ($expmtType =~ /two-colour array/) { $color = "2colour_" ;} 

my $RNA = "mrna";
if (exists $H_miRnaList{$expAcc}) { $RNA = "microrna" ;}

my $analysisType ;
if ($differential) { $analysisType = "differential" ;} 
if ($baseline) { $analysisType = "baseline" ;}

## At the moment, we exclude microarray baseline experiments, so die if we have one of those
if ( $type eq "microarray" && $baseline) { 
	  &log($logFileHandle, "[ERROR] $expAcc - $type experiment cannot be analysed as baseline. Refused for now.") ; die ; }

## For each experiment type, make sure I've got the appropriate bits
if ($type =~ /rnaseq/ && $RNA ne '' && $analysisType ne '') {
	$experimentType = "${type}_${color}${RNA}_${analysisType}" ;
    if ($debug) { &log($logFileHandle, "Exepriment (Atlas) type is $experimentType (from \"$expmtType\")") ; }
} elsif ($type =~ /microarray/ && $RNA ne '' && $color ne '') {
    $experimentType = "${type}_${color}${RNA}_${analysisType}" ;
    if ($debug) { &log($logFileHandle, "Exepriment (Atlas) type is $experimentType (from \"$expmtType\")") ; }
} else { &log($logFileHandle, "[ERROR] $expAcc - Cannot get the experiment type: type:$type color:$color RNA:$RNA [from $expmtType]\n") ; die ; } 



############################
## Factor Value parsing
## For each organism and each array design (usually: 1 of each only),
## Parse the factor values: 
#	- >= 3 ior 2 replicates	  -- delete the Factor Value if not true 
#	- reference		          -- from the list of references generated from the config file, or passed as argument
#	- forbidden Factor Values -- from the kill list generated from the config file, or passed as argument; delete Factor Value if true
my $singleFactor ;
my $configurationTag = 0 ;
$noReferenceError = "Candidate reference values for $factorType: ";
if ($debug) { &log($logFileHandle, "[DEBUG] Parsing values collected in Magetab module\n") ; }

## Number of arrays
#  - RNA_seq experiment: always one
#  - microarray experiment: possibly more than one
my $arrayNumber = scalar keys %H_eFactorValues2runIDs ;
if ($debug) { &log($logFileHandle, "[DEBUG] Number of arrays in that experiment (initially): $arrayNumber") ; }

my %H_referenceArray ; #store reference for each couple array/organism

if ($debug) { &log($logFileHandle, "\n[DEBUG] ===== PARSING the DATA (cleaning step) =====\n") ; }
foreach my $array (sort keys %H_eFactorValues2runIDs) {
    if ($debug) { &log($logFileHandle, "[DEBUG] Array is $array ($factorType)") ; }

	## Report error when testing the replicates, reference etc. 
	my $warningMessage = "" ;

	foreach my $organism (sort keys %{$H_eFactorValues2runIDs{$array}}) {
		if ($debug) { &log($logFileHandle, "[DEBUG]\tSpecies is $organism\n") ; }

		foreach my $time (sort {$a<=>$b} keys %{$H_eFactorValues2runIDs{$array}{$organism}}) {
			if ($debug) { &log($logFileHandle, "[DEBUG]\tTime is $time\n") ; }

			my $arrayOrgTime = "$array;$organism;$time" ;

			foreach my $FVlist (keys %{$H_eFactorValues2runIDs{$array}{$organism}{$time}}) {
	
				## Differential analysis only
				if ($differential) {

					## Split FVlist, as could be a string of several FVs (X-factor experiments) 
					my @A_FV = split ("\t\t", $FVlist) ;

					for my $FV (@A_FV) {

						## Remove forbidden factor value (e.g. 'individual')
						if (exists $H_config{"FACTOR_VALUE_KILL"}{lc($FV)}) { delete $H_eFactorValues2runIDs{$array}{$organism}{$time}{$FVlist} ; next ; } 		
						$noReferenceError .= " '$FV' ";

						## Test for reference - if already one, die loudly
						# (case insensitive: lc only)
						if (exists $H_config{"REFERENCE"}{lc($FV)}) {
							
							## If only 1 factor value: we only want 1 reference
							if (scalar(@A_FV) == 1) {

								## Set a single factor flag
								$singleFactor = 1 ; 

								if (!exists $H_referenceArray{$arrayOrgTime}) {
									$H_referenceArray{$arrayOrgTime}{$FVlist} = 1 ;
									if ($debug) { &log($logFileHandle, "[DEBUG] Adding reference $FVlist (with $arrayOrgTime) to the list [1]\n") ; }
								} else { 
									$warningMessage .= "More than one reference: $H_referenceArray{$arrayOrgTime} and $FV. " ;
								}
							}

							## If more than 1 factor value: we accept multiple references
							else {
								$H_referenceArray{$arrayOrgTime}{$FVlist} = 1 ;
								if ($debug) { &log($logFileHandle, "[DEBUG] Adding reference $FVlist (with $arrayOrgTime) to the list [2]\n") ; } 
							}
						}
					}
				}

				## Test for statistical replicates
				# Differential: at least 3
				# Baseline: at least 2, unless -noreplicate tag (for specific experiments that we really want in Atlas)
				#
				# If less replicates than required, then delete entries for this F.V. in:
				#    - %H_eFactorValues2runIDs{$array}{$organism}{$time}{$FV}
				#    - %H_referenceArray{$array;$organism;$time}{$FV}
				#    - %H_assayGroups1Difference{$array;$organism;$time}{$FV}  and  %H_assayGroups1Difference{$array;$organism;$time}{x}{$FV}
				my $replicateCount = scalar @{$H_eFactorValues2runIDs{$array}{$organism}{$time}{$FVlist}} ;
				if ($debug) { &log($logFileHandle, "[DEBUG] Replicate number: $replicateCount\n") ; }

				my $deleteFlag = 0 ;
				if ($differential  && ($replicateCount < 3)) { $deleteFlag = 1 ; }
				if ($baseline && ($replicateCount < 2) && !$noreplicate) { $deleteFlag = 1 ; }
				if ($deleteFlag eq 1) {
					delete $H_eFactorValues2runIDs{$array}{$organism}{$time}{$FVlist} ; 	

					#'if exists' required otherwise autovivification with the 'delete'
					if ($H_referenceArray{$arrayOrgTime} && exists $H_referenceArray{$arrayOrgTime}{$FVlist}) { delete $H_referenceArray{$arrayOrgTime}{$FVlist} ; } 
					foreach my $assay (keys %{$H_assayFactorValues{$arrayOrgTime}}) { 
						if ($H_assayFactorValues{$arrayOrgTime}{$assay} eq $FVlist) { 
							&log($logFileHandle, "[DEBUG] Deleting H_assayFactorValues{$arrayOrgTime}{$assay}: $H_assayFactorValues{$arrayOrgTime}{$assay}\n") ;
							delete $H_assayFactorValues{$arrayOrgTime}{$assay} ;
						}
					}
					if ($debug) { &log($logFileHandle, "[DEBUG] Delete H_eFactorValues2runIDs{$array}{$organism}{$time}{$FVlist} due to lack of replicates\n") ; }
					if ($debug) { &log($logFileHandle, "[DEBUG] Delete H_referenceArray{$arrayOrgTime}{$FVlist} due to lack of replicates\n") ; }
					if ($debug) { &log($logFileHandle, "[DEBUG] Delete H_assayFactorValues{$arrayOrgTime}{x} = $FVlist due to lack of replicates\n") ; }
				}	
			}

			## For differential only 
			# Test if reference Factor Value found 
			if ($differential && $singleFactor) {

				if (!exists $H_referenceArray{$arrayOrgTime}) {
					$warningMessage .= "No reference: $noReferenceError. Add it to $conf or add it temporarily with -ref \"new reference word\"." ;
				}
			
				## Check number of Factor Value left (with >= 3 replicates)
				# ... need at least two sets of Factor Values
				# ... with maximum 1 factor different (e.g. in case of X-factors)

				my %H_storeDiff ; #record the number of differences between 2 factor values 
				if (keys %{$H_eFactorValues2runIDs{$array}{$organism}{$time}} < 2) {
					$warningMessage .= "Less than 2 factor values with at least 3 replicates. " ;
					
					## If time factor values, this is a situation we want to record, 
					# As it might influence what we do later
					# Record the time and the whole set of values
					if ($time ne "NOTIME") {
						if ($debug) { &log($logFileHandle, "Time factor ($time) and <2 F.V. set with 3 reps. Archive for later!\n") ; } 
						$H_timePointsNotEnoughFvSets{$array}{$organism}{$time} = 1 ;

						## Single F.V. for that array;organism;time
						# But need to store a variable, as opposed to a single value hash
						# (easier to retrieve)
						foreach my $fv (keys %{$H_eFactorValues2runIDs{$array}{$organism}{$time}}) {
							$H_timePointsNotEnoughFvSetsFV2runIDs{$array}{$organism}{$time} = $fv ;
						}
					}	
				} else { 
					foreach my $FVtest1 (keys %{$H_eFactorValues2runIDs{$array}{$organism}{$time}}) {
						foreach my $FVtest2 (keys %{$H_eFactorValues2runIDs{$array}{$organism}{$time}}) {
							$H_storeDiff{compareLists($FVtest1,$FVtest2)}++ ;
						}
					}

					if (!exists $H_storeDiff{1}) {
                		$warningMessage .= "No factor value pair with only 1 difference. " ;
					}

					if ($H_storeDiff{1} < 2) {
						$warningMessage .= "Less than 2 factor value pairs with only 1 difference. " ;
					}
				}
			}


			## Cannot generate contrast file
			# Warn, don't die
			if ($warningMessage ne "") {

				## If time points, there are some errors we accept: 
				#   - no reference
				#   - no factor value pair with only 1 difference
				#   - less than 2 factor values with at least 3 replicates
				#   - less than 2 factor value pairs with only 1 difference
				#
				# And some that are deadly: 
				#   - more than one reference (per array, per org., per time point)

				## So only delete if NOTIME, or if 'more than one reference'
				if ( ($time eq "NOTIME") || ($warningMessage =~ /More than one reference/)) {
                	&log($logFileHandle, "\n[WARNING] XML configuration file cannot be generated for $expAcc :: $organism :: $time :: $array: $warningMessage\n") ;
					if ($debug) { &log($logFileHandle, "[DEBUG] Delete $array ; $organism ; $time ...\n") ; }
					delete $H_eFactorValues2runIDs{$array}{$organism}{$time} ; 
				}
			}
		}
	}
}
if ($debug) { &log($logFileHandle, "[DEBUG] Finish reading Magetab files\n") ; }



############################
### Now go through the hash again and generate assay_groups
## This can only be done once I have a whole view of what's left
## after removal of factor values with not enough replicates, or reference.
#
## Generate assay_groupS
# Get the reference FV: 
#   - other FVs that differs with 1 factor?
#   --- if no: no contrast
#   --- if yes, build the assayGroups from there
#
# If no reference: 
#  - do a all vs. all.
#  --- so build an assayGroup with all of them
#
# Data stored in: 
#  - references: %H_referenceArray{array;organism;time}{FVlist}
#  - %H_assayFactorValues{array;organism;time}{assay} = list F.V.
#  - ... and compare those groups
#######
if ($debug) { &log($logFileHandle, "\n[DEBUG] ===== PARSING the DATA (<assay_groups>) =====\n") ; }
undef %H_assayGroups1Difference ;
foreach my $array (sort keys %H_eFactorValues2runIDs) {
	
	foreach my $organism (sort keys %{$H_eFactorValues2runIDs{$array}}) {
		if ($debug) { &log($logFileHandle, "[DEBUG]\tSpecies is $organism\n") ; } 
		
		foreach my $time (sort keys %{$H_eFactorValues2runIDs{$array}{$organism}}) {
			if ($debug) { &log($logFileHandle, "[DEBUG]\tTime is $time\n") ; }

			my $arrOrgTime = "$array;$organism;$time" ;

			## Compare lists of F.V.
			#  ... store them both ways $H_assayGroups1Difference{A}{B} & {B}{A}
			#  .... easier to check for existence later on.
			foreach my $assayName1 (sort keys %{$H_assayFactorValues{$arrOrgTime}}) {
				foreach my $assayName2 (sort keys %{$H_assayFactorValues{$arrOrgTime}}) {
					if ($assayName1 ne $assayName2) {
						my $FVlist1 = $H_assayFactorValues{$arrOrgTime}{$assayName1} ;
						my $FVlist2 = $H_assayFactorValues{$arrOrgTime}{$assayName2} ;

						## Compare those 2 lists, 
						# If 1 difference: same assay_group
						if (compareLists($FVlist1,$FVlist2) == 1) {
							$H_assayGroups1Difference{$arrOrgTime}{$FVlist1}{$FVlist2} = 1 ;
						}
					}
				}
			}

			## For differential
			my $groupCounter = &valuesNumberHashOfHash(\%H_inAssayGroup) +1 ;
			if ($debug) { &log($logFileHandle, "Group counter, $array, starting at $groupCounter\n") ; }

			if ($differential) {
				## Reference: ref. vs. ALL
				if (exists $H_referenceArray{$arrOrgTime}) {
					$differentialTag = "REFERENCE" ;
						
					if ($debug) { &log($logFileHandle, "=-=-=-=-= Making the groups (ref.)...\n") ; }
					foreach my $referenceFV (sort keys %{$H_referenceArray{$arrOrgTime}}) {
						if ($debug) { &log($logFileHandle, "[DEBUG] And reference is: $referenceFV\n") ; }

						## If present in %H_assayGroups1Difference
						# i.e. linked to another assay_group with exactly 1 difference in F.V.
						if (exists $H_assayGroups1Difference{$arrOrgTime}{$referenceFV}) {

							## If not already in an assay_group,
							# assign the reference to one
							if (!exists $H_inAssayGroup{$arrOrgTime}{$referenceFV}) {
								$H_inAssayGroup{$arrOrgTime}{$referenceFV} = $groupCounter ;
								$groupCounter++ ;
								if ($debug) { &log($logFileHandle, "[DEBUG] [Assay_group $H_inAssayGroup{$arrOrgTime}{$referenceFV}] $referenceFV\n") ; }
							}
				
							foreach my $FV (sort keys %{$H_assayGroups1Difference{$arrOrgTime}{$referenceFV}}) {

								## If not already in an assay_group,
								# assign the F.V. to one
								if (!exists $H_inAssayGroup{$arrOrgTime}{$FV}) {
									$H_inAssayGroup{$arrOrgTime}{$FV} = $groupCounter ;
									$groupCounter++ ;
									if ($debug) { &log($logFileHandle, "[DEBUG] [Assay_group $H_inAssayGroup{$arrOrgTime}{$FV}] $FV\n") ; }
								}

								## Mark this F.V. pair as seen
								$H_assayGroups1Difference{$arrOrgTime}{$referenceFV}{$FV} = 0 ; #1 if pair never seen    
						
								## And delete the opposite pair
								## Avoid parsing the same pair twice as A,B and B,A
								delete $H_assayGroups1Difference{$arrOrgTime}{$FV}{$referenceFV} ; 
							}
						}
					}

				## If no reference: ALL vs. ALL
				# assayGroups order is at random
				} else {
					if ($debug) { &log($logFileHandle, "=-=-=-=-= Making the groups (no ref.)...\n") ; }
					$differentialTag = "NOREF" ;
					foreach my $FV1 (sort keys %{$H_assayGroups1Difference{$arrOrgTime}}) {
	
						## If not already in an assay_group,
						# assign the $FV1 to one
						if (!exists $H_inAssayGroup{$arrOrgTime}{$FV1}) {
							$H_inAssayGroup{$arrOrgTime}{$FV1} = $groupCounter ;
							$groupCounter++ ;
							if ($debug) { &log($logFileHandle, "[DEBUG] [Assay_group $H_inAssayGroup{$arrOrgTime}{$FV1}] $arrOrgTime $FV1\n") ; }
						}

						foreach my $FV2 (sort keys %{$H_assayGroups1Difference{$arrOrgTime}{$FV1}}) {

							## If not already in an assay_group,
							# assign the $FV2 to one
							if (!exists $H_inAssayGroup{$arrOrgTime}{$FV2}) {
								$H_inAssayGroup{$arrOrgTime}{$FV2} = $groupCounter ;
								$groupCounter++ ;
								if ($debug) { &log($logFileHandle, "[DEBUG] [Assay_group $H_inAssayGroup{$arrOrgTime}{$FV2}] $arrOrgTime $FV2\n") ; }
							}

							## Mark this F.V. pair as seen
							if ($debug) { &log($logFileHandle, "[DEBUG] Comparing 1 diff: {$arrOrgTime}{$FV1}{$FV2}\n") ; } 
							$H_assayGroups1Difference{$arrOrgTime}{$FV1}{$FV2} = 0 ; #1 if pair never seen   

							## And delete the opposite pair
							## Avoid parsing the same pair twice as A,B and B,A
							if ($debug) { &log($logFileHandle, "[DEBUG] Deleting oppposite pair: {$arrOrgTime}{$FV2}{$FV1}\n") ; }
							delete $H_assayGroups1Difference{$arrOrgTime}{$FV2}{$FV1} ;
						}
					}
				}
			}

			## For baseline
			## No reference, so assayGroups order is at random
			if ($baseline) {
				foreach my $FV (sort keys %{$H_eFactorValues2runIDs{$array}{$organism}{$time}}) {
					$H_inAssayGroup{$arrOrgTime}{$FV} = $groupCounter ;
					$groupCounter++ ;
				}
			}
		}
	}
}

## If debug:
# Print assay_groups generated - no reference issue for now.
if ($debug) {
	&log($logFileHandle, "=-=-=-=-= Possible comparisons (no ref., 2-sided)...\n") ;
	foreach my $aOT (sort keys %H_assayGroups1Difference) {
		foreach my $list1 (sort keys %{$H_assayGroups1Difference{$aOT}}) {
			foreach my $list2 (sort keys %{$H_assayGroups1Difference{$aOT}{$list1}}) {
				&log($logFileHandle, "[$aOT] $list1 <-> $list2 ($H_assayGroups1Difference{$aOT}{$list1}{$list2})\n") ;
			} &log($logFileHandle, "--\n") ;
		}
	}
}

## If factor time and no time point with >=2 sets of Factor Values
# ... compare time point to each other (if 1 factor value difference)
#
## For differential only, 
# %H_timePointsNotEnoughFvSets has been created in a differential context 
# so no need to test for it. 
if (keys %H_timePointsNotEnoughFvSets > 0) {

	if ($debug) { &log($logFileHandle, "[DEBUG] Do a time series comparison...\n") ; }

	## Creating an array with an ordered list of time points
	# At the moment I'm assuming they are in the same unit and I can just sort numerically
	foreach my $array (sort keys %H_timePointsNotEnoughFvSets) {
    	foreach my $organism (sort keys %{$H_timePointsNotEnoughFvSets{$array}}) {
			my $cpt = 0 ;

			## Generate an array with time in order
			foreach my $time (sort {$a<=>$b} keys %{$H_timePointsNotEnoughFvSets{$array}{$organism}}) {
				$A_timePointsNotEnoughFvSets[$cpt] = $time ;
				$cpt++ ;
			}

			## Generate <assay_group>
			# Smallest time point: reference or not?
			my $firstTimePoint = $A_timePointsNotEnoughFvSets[0] ;
			my $firstFactorValue = $H_timePointsNotEnoughFvSetsFV2runIDs{$array}{$organism}{$firstTimePoint} ;
			my $arrOrgTime = "$array;$organism;$firstTimePoint" ;
			my $groupCounter = 1 ;
			if ($debug) { &log($logFileHandle, "[DEBUG] First time point: $firstTimePoint and associated f.v. $firstFactorValue\n") ; }

           	## If reference: 
			#  - compare T0 to T1, T2 etc.
			#  - IF 1 diff. between f.v. sets
			#  - create <assay_group>
			# Data in $H_timePointsNotEnoughFvSetsFV2runIDs{$array}{$organism}{$time}
			if (exists $H_referenceArray{$arrOrgTime}{$firstFactorValue}) {
				
				for my $a (1..$#A_timePointsNotEnoughFvSets) { #[0] is the reference
					my $timePoint = $A_timePointsNotEnoughFvSets[$a] ;
				
					# Factor value at that time point
					my $FV = $H_timePointsNotEnoughFvSetsFV2runIDs{$array}{$organism}{$timePoint} ;
					if ($debug) { &log($logFileHandle, "[DEBUG] Testing with ref. - FV is $FV and time is $timePoint\n") ; }

					## One difference between FVref and FVtimeT?
					#  If so, create an assay_group
					## %H_assayGroups1Difference{$arrOrgTime}{$FVa}{$FVb} doesn't contain F.V. pair ACROSS time points
					#  So do the comparison on the fly
                    ## Since doing a time series comparison, 
					#  FV in $H_inAssayGroup{$arrOrgTime}{$FV} needs to contains time 
					#  e.g., we can have t0 vs. t1, ref. vs. FV-a
					#               and t0 vs. t2, ref. vs. same FV-a
					# ... which means it'll need to be removed when printing 
					if (compareLists($firstFactorValue,$FV) == 1) {
						if (!exists $H_inAssayGroup{$arrOrgTime}{$firstFactorValue}) { 
							$H_inAssayGroup{$arrOrgTime}{$firstFactorValue} = $groupCounter ;
							$groupCounter++ ;
							if ($debug) { &log($logFileHandle, "[DEBUG] Storing assay_group REF H_inAssayGroup{$arrOrgTime}{$firstFactorValue} ; assay_group $H_inAssayGroup{$arrOrgTime}{$firstFactorValue}\n") ; }
						}

						$FV = $FV.";;;".$timePoint ;

						$H_inAssayGroup{$arrOrgTime}{$FV} = $groupCounter ;
						$groupCounter++ ;
						if ($debug) { &log($logFileHandle, "[DEBUG] Storing assay_group H_inAssayGroup{$arrOrgTime}{$FV} ; assay_group $H_inAssayGroup{$arrOrgTime}{$FV}\n") ; }

						# Populating a %H similar to %H_assayGroups1Difference
						# ... but specific for time series
						$H_assayGroups1DifferenceTimeSeries{$arrOrgTime}{$firstFactorValue}{$FV} = 1 ;
					}
				}
			}

			## No reference: 
			#   - compare T0 to T1, T1 to T2, T2 to T3, etc.
			#   - IF 1 diff. between f.v. sets
			#   - create <assay_group>
			else { 
				if ($debug) { &log($logFileHandle, "[DEBUG] NO REFERENCE at lower time point! Doing t1 vs. t2 ....\n") ; }

				for my $a (0..$#A_timePointsNotEnoughFvSets-1) {

					# Get the time points
					my $timePoint = $A_timePointsNotEnoughFvSets[$a] ;
					my $nextTimePoint = $A_timePointsNotEnoughFvSets[$a+1] ;

					# Get the cognate F.V.
					my $FV = $H_timePointsNotEnoughFvSetsFV2runIDs{$array}{$organism}{$timePoint} ;
					my $nextFV = $H_timePointsNotEnoughFvSetsFV2runIDs{$array}{$organism}{$nextTimePoint} ;

					if ($debug) { &log($logFileHandle, "[DEBUG] Time points are: $timePoint [$a] & $nextTimePoint [$a+1]. F.V. are $FV & $nextFV\n") ; } 

					## One difference between F.V. sets?
					#  If so, create an assay_group
					## %H_assayGroups1Difference{$arrOrgTime}{$FVa}{$FVb} doesn't contain F.V. pair ACROSS time points
					#  So do the comparison on the fly
					## Since doing a time series comparison, 
					#  FV in $H_inAssayGroup{$arrOrgTime}{$FV} needs to contains time 
					#  e.g., we can have t0 vs. t1, ref. vs. FV-a 
					#                and t0 vs. t2, ref. vs. same FV-a
					#  ... which means it (timw) will need removing when printing              
					if (compareLists($FV, $nextFV) == 1) {

						if (!exists $H_inAssayGroup{$arrOrgTime}{$FV}) {
							$FV = $FV.";;;".$nextTimePoint ;
							$H_inAssayGroup{$arrOrgTime}{$FV} = $groupCounter ;
							$groupCounter++ ;
						}
						if (!exists $H_inAssayGroup{$arrOrgTime}{$nextFV}) {
							$H_inAssayGroup{$arrOrgTime}{$nextFV} = $groupCounter ;
							$groupCounter++ ;
						}

					   ## Populating a %H similar to %H_assayGroups1Difference
					   # ... but specific for time series
					   $H_assayGroups1DifferenceTimeSeries{$arrOrgTime}{$FV}{$nextFV} = 1 ;
					}
				}
			}
		}
	}
}


############################
#### Now go through the hash again, and print
## In some cases (e.g. : Time factor value, multi-organism), this can only be done 
## once I have all the <assay_groupS> and all the <contrasts>

## Array number
## 1. In the hash, delete any array with no organism left
foreach my $arrCheck (keys %H_eFactorValues2runIDs) {

	foreach my $orgCheck (keys %{$H_eFactorValues2runIDs{$arrCheck}}) {
		## If no F.V.: delete organism
		if (scalar keys %{$H_eFactorValues2runIDs{$arrCheck}{$orgCheck}} == 0) {
			delete $H_eFactorValues2runIDs{$arrCheck}{$orgCheck} ;
		}
	}   
												    
	## If no organism: delete array
	if (scalar keys %{$H_eFactorValues2runIDs{$arrCheck}} == 0) {
		delete $H_eFactorValues2runIDs{$arrCheck} ;
	}
}

## 2. Calculate the array number 
# Required later to decide what to print
my $arrayNumber = scalar keys %H_eFactorValues2runIDs ;
if ($debug) { &log($logFileHandle, "[DEBUG] Number of arrays in that experiment (after parsing): $arrayNumber\n") ; }

## Print the XML
if ($debug) { &log($logFileHandle, "\n[DEBUG] ===== PRINTING the XML =====\n") ; }

## Is there anything to print?
if ($debug) {
	my $anyAssayGrp = keys %H_inAssayGroup ;
	&log($logFileHandle, "[DEBUG] Assay group left for printing: $anyAssayGrp\n") ;
}


if ( (keys %H_inAssayGroup) != 0 ) {

	## Open file to write to
	$XML = IO::File->new(">$outfileXML");

	## Use the XML::Writer module to print
	$writer = XML::Writer->new(OUTPUT => $XML, DATA_MODE => 1, DATA_INDENT => 4);

	## Begin configuration XML, add experiment type.
	$writer->startTag("configuration", "experimentType" => $experimentType);

	## For each array	
	foreach my $array (sort keys %H_eFactorValues2runIDs) {

		$writer->startTag("analytics");

		## Array_design element - only if experiment has an array
		if ($array ne "0") { $writer->dataElement("array_design" => $array) ; }

		## <assay_groupS> element
		$writer->startTag("assay_groups");
		if ($debug) { &log($logFileHandle, "[DEBUG] Opening <assay_groupS>\n") ; }

		foreach my $arrOrgTime (sort keys %H_inAssayGroup) {
			if ($arrOrgTime =~ /$array/) {

				## Individual <assay_group> elements
				foreach my $factVal (sort { $H_inAssayGroup{$arrOrgTime}{$a} <=> $H_inAssayGroup{$arrOrgTime}{$b} } keys %{$H_inAssayGroup{$arrOrgTime}}) {
					if ($debug) { &log($logFileHandle, "[DEBUG] Array, organism and time: $arrOrgTime\n") ; } 
					if ($debug) { &log($logFileHandle, "[DEBUG] [$H_inAssayGroup{$arrOrgTime}{$factVal}] F.V. is $factVal\n") ; }

					## Label
					my $label = "" ;

					## If time series, 
					#  ... then information is stored differently in %H_inAssayGroup
					#  ... so get the appropriate $arrOrgTime and $factVal (->$arrOrgTime_2 & $factVal_2)
					#  ... to query %H_factorValue2factorType 
					my $arrOrgTime_2 = $arrOrgTime ;
					my $factVal_2 = $factVal ;

					## If time series
					if ($factVal =~ /^(.+?);;;(.+?)$/) { 	
						$factVal_2 = $1 ;
						my $newTime = $2 ;
						my @A_arrOrgTime = split (";",$arrOrgTime) ;
						$arrOrgTime_2 = "$A_arrOrgTime[0];$A_arrOrgTime[1];$newTime" ;
						if ($debug) { &log($logFileHandle, "Time series, new value for arrOrgTime_2 ($arrOrgTime_2) and factVal_2 ($factVal_2)\n") ; } 
					}

					foreach my $ft (sort keys %{$H_factorValue2factorType{$arrOrgTime_2}{$factVal_2}}) {
						##$label .= "$ft: $H_factorValue2factorType{$arrOrgTime}{$factVal}{$ft}; " ;  #label="genotype: stwl bam double mutant"
						$label .= "$H_factorValue2factorType{$arrOrgTime_2}{$factVal_2}{$ft}; " ;  #label="stwl bam double mutant"
					}

					## Remove the trailing
					$label =~ s/; $// ;

					## Write element
					$writer->startTag("assay_group", "id" => "g$H_inAssayGroup{$arrOrgTime}{$factVal}", "label" => $label) ;
					if ($debug) { &log($logFileHandle, "[DEBUG] Opening <assay_group>\n") ; }
					if ($debug) { &log($logFileHandle, "assay_group: $factVal - $label\n") ; }

					my ($array,$organism,$time) = split (";", $arrOrgTime_2) ;	
					foreach my $names (@{$H_eFactorValues2runIDs{$array}{$organism}{$time}{$factVal_2}}) {
						if ($debug) { &log($logFileHandle, "\t- assay name: $names\n") ; }

						if ($names ne '') {	
							## Technical replicates, if any
							my $tech_rep ;
							if ($H_TechRepsGroup{$names}) { $writer->dataElement("assay" => $names, "technical_replicateid" => $H_TechRepsGroup{$names}) ; }
							else { $writer->dataElement("assay" => $names) ; }
						}
					}
					$writer->endTag("assay_group") ;
					if ($debug) { &log($logFileHandle, "[DEBUG] Closing <assay_group>\n") ; }
				}
			}
		}

		## Close <assay_groupS> element
		$writer->endTag("assay_groups");
		if ($debug) { &log($logFileHandle, "[DEBUG] Closing <assay_groupS>\n") ; }

		## Contrast element (differential only)
		# There are 2 sections:
		#    - normal contrasts, most common
		#    (ref. vs. rest, or ALL vs. ALL)
		#    - time series
		#    (T0 vs. the rest, or Tn vs. Tn+1)
		#
		## Classic contrasts (NOT for time series)
		# Reference can have any ID but: 
		#   - in the contrast ID the reference ID comes first
		#   - in the contrast name the reference name comes last
		if ($differential && keys %H_assayGroups1DifferenceTimeSeries == 0) {

			## Open the general <contrasts> element
			$writer->startTag("contrasts") ;
			if ($debug) { &log($logFileHandle, "[DEBUG] Opening <contrastS> in classic contrast mode\n") ; }
		
			## Contrasts with reference(s):
			#    - for each reference in %H_referenceArray,
			#    - get the paired F.V. from %H_assayGroups1Difference
			if ($singleFactor || $differentialTag eq "REFERENCE") {
			            
				if ($debug) { &log($logFileHandle, "=-=-=-=-= Print the contrast\n") ; }

				if ($debug) { &log($logFileHandle, "Going through the references...\n") ; }
				foreach my $arrOrgTime (sort %H_referenceArray) {
					foreach my $refFV (sort keys %{$H_referenceArray{$arrOrgTime}}) {
						if ($arrOrgTime =~ /$array/) {
							if ($debug) { &log($logFileHandle, "[DEBUG] Ref. is $refFV\n") ; }
							my $refGroup = $H_inAssayGroup{$arrOrgTime}{$refFV} ;

							foreach my $nonRefFV (sort keys %{$H_assayGroups1Difference{$arrOrgTime}{$refFV}}) {
								my $nonRefGroup = $H_inAssayGroup{$arrOrgTime}{$nonRefFV} ;

								if ($debug) { &log($logFileHandle, "REFERENCE: $refFV [$refGroup] vs. $nonRefFV [$nonRefGroup]\n") ; }
								$writer->startTag("contrast", "id" => "g${refGroup}_g$nonRefGroup") ;

								if ($debug) { &log($logFileHandle, "[DEBUG] Printing <contrast> - reference\n") ; }

								## Replace the \t\t separators by ;
								## So match the label in <assay_group id="gZ" label="xxx; yyy"> 
								$nonRefFV =~ s/\t\t/; /g ;
								$refFV =~ s/\t\t/; /g ;

								## If F.V. time, add time to the contrast name
								my $contrastTime = "" ;
								my ($arr, $org, $time) = split(";", $arrOrgTime) ;  
								if ($time ne "NOTIME") { $contrastTime = " at $time" ; }

								if ($arrayNumber == 1) { $writer->dataElement("name" => "'$nonRefFV' vs '$refFV'$contrastTime") ; }

								## If > 1 array
								## user friendly name for the array
								else {
									if ($H_arrayIDs2arrayNames{$array} ne '') {	
										$writer->dataElement("name" => "'$nonRefFV' vs '$refFV'$contrastTime on '$H_arrayIDs2arrayNames{$array}'") ;
                    				} else { &log($logFileHandle, "[ERROR] $expAcc - No user-friendly name for array $array") ; die ; }
								}
								$writer->dataElement("reference_assay_group" => "g$refGroup") ;
								$writer->dataElement("test_assay_group" => "g$nonRefGroup") ; 
								$writer->endTag("contrast") ;
							}
						}
					}
				}

			## If no reference: all vs. all
			# We need to choose a random reference:
			# 	- $FV1 will be the reference, and $FV2 will be the non-reference
			# 	- same order as for reference contrast:
			# 		"contrast", "id" => "g-refgroup_g-nonrefgroup"
			# 	but: 
			#		"name" => "'non-ref' vs 'ref'" 
			} else {
				foreach my $arrOrgTime (sort keys %H_assayGroups1Difference) {
					foreach my $FV1 (sort keys %{$H_assayGroups1Difference{$arrOrgTime}}) {

						## Get FV2 assay_group
						my $fv1Group = $H_inAssayGroup{$arrOrgTime}{$FV1} ;
						foreach my $FV2 (sort keys %{$H_assayGroups1Difference{$arrOrgTime}{$FV1}}) {

							if ($arrOrgTime =~ /$array/) {

								if ($debug) { &log($logFileHandle, "Testing $FV1 & $FV2 at $arrOrgTime\n") ; }

								## If pair never seen
								if (!exists $H_assayGroupDone{$arrOrgTime}{$FV1}{$FV2}) {

        	   	            		## Store that pair so won't do it again (A,B or B,A)
									$H_assayGroupDone{$arrOrgTime}{$FV1}{$FV2} = 1 ;
									$H_assayGroupDone{$arrOrgTime}{$FV2}{$FV1} = 1 ;

									## Get FV2 assay_group
									my $fv2Group = $H_inAssayGroup{$arrOrgTime}{$FV2} ;
 
									$writer->startTag("contrast", "id" => "g${fv1Group}_g$fv2Group") ;
									if ($debug) { &log($logFileHandle, "[DEBUG] Printing <contrast> - $FV2 (g$fv2Group) vs. $FV1 (g${fv1Group}) - no ref.\n") ; }

									## If 1 array, or RNA_seq
									# Replace the \t\t separators by ;
									# So match the label in <assay_group id="gZ" label="xxx; yyy"> 
									$FV1 =~ s/\t\t/; /g ;
									$FV2 =~ s/\t\t/; /g ;

									## If F.V. time, add time to the contrast name
									my $contrastTime = "" ;
									my ($arr, $org, $time) = split(";", $arrOrgTime) ; 	
									if ($time ne "NOTIME") { $contrastTime = " at $time" ; }
										
									if ($arrayNumber == 1) { $writer->dataElement("name" => "'$FV2' vs '$FV1'$contrastTime") ; }
						
									## If > 1 array,
									# user friendly name for the array
									else {
										if ($H_arrayIDs2arrayNames{$array} ne '') {
											$writer->dataElement("name" => "'$FV2' vs '$FV1'$contrastTime on '$H_arrayIDs2arrayNames{$array}'") ;
                                    	} else { &log($logFileHandle, "[ERROR] $expAcc - No user-friendly name for array $array") ; die ; }
									}
									$writer->dataElement("reference_assay_group" => "g$fv1Group") ;
									$writer->dataElement("test_assay_group" => "g$fv2Group") ;
									$writer->endTag("contrast") ;
								}
							}
						}
					}
				}
			}
			if ($debug) { &log($logFileHandle, "[DEBUG] Closing <contrastS>\n") ;  }
			$writer->endTag("contrasts") ;
		}


		## Contrast for time series
		# Reference can have any ID but: 
		#   - in the contrast ID the reference ID comes first
		#   - in the contrast name the reference name comes last
		if (keys %H_assayGroups1DifferenceTimeSeries > 0) {

			## Open the general <contrasts> element
			$writer->startTag("contrasts") ;
			if ($debug) { &log($logFileHandle, "[DEBUG] Opening <contrastS> in time series mode\n") ; }

			## Contrasts with reference(s):
			#   - for each reference in %H_referenceArray,
			#   - get the paired F.V. from %H_assayGroups1Difference
			if ($differentialTag eq "REFERENCE") {
				if ($debug) { &log($logFileHandle, "=-=-=-=-= Print the contrast (time series)\n") ; }
				if ($debug) { &log($logFileHandle, "Going through the references...\n") ; }
			
				foreach my $arrOrgTime (sort %H_referenceArray) {
					foreach my $refFV (sort keys %{$H_referenceArray{$arrOrgTime}}) {
						if ($arrOrgTime =~ /$array/) {
							if ($debug) { &log($logFileHandle, "[DEBUG] Ref. is $refFV\n") ; }

							foreach my $nonRefFV (sort keys %{$H_assayGroups1DifferenceTimeSeries{$arrOrgTime}{$refFV}}) {
								my $refGroup = $H_inAssayGroup{$arrOrgTime}{$refFV} ;
								my $nonRefGroup = $H_inAssayGroup{$arrOrgTime}{$nonRefFV} ;
							
								if ($debug) { &log($logFileHandle, "REFERENCE: $refFV [$refGroup] vs. $nonRefFV [$nonRefGroup]\n") ; }								
								$writer->startTag("contrast", "id" => "g${refGroup}_g$nonRefGroup") ;
								
								if ($debug) { &log($logFileHandle, "[DEBUG] Printing <contrast> - reference\n") ; }

								## Remove the ;;;
								$nonRefFV =~ s/;;;/; /g ;
								$refFV =~ s/;;;/; /g ;

								## Replace the \t\t separators by ;
								# ... so match the label in <assay_group id="gZ" label="xxx; yyy"> 
								$nonRefFV =~ s/\t\t/; /g ;
								$refFV =~ s/\t\t/; /g ;
								if ($arrayNumber == 1) { $writer->dataElement("name" => "'$nonRefFV' vs '$refFV'") ; }

								## If > 1 array,
								# user friendly name for the array
								else {
									if ($H_arrayIDs2arrayNames{$array} ne '') {
										$writer->dataElement("name" => "'$nonRefFV' vs '$refFV' on '$H_arrayIDs2arrayNames{$array}'") ;
									} else { die "[ERROR] $expAcc - No user-friendly name for array $array\n" ; }
								}
                                $writer->dataElement("reference_assay_group" => "g$refGroup") ;
								$writer->dataElement("test_assay_group" => "g$nonRefGroup") ;
								$writer->endTag("contrast") ;
							}
						}
					}
				}

			## If no reference: all vs. all
			## We need to choose a random reference:
			#    - $FV1 will be the reference, and $FV2 will be the non-reference
			#    - same order as for reference contrast:
			#        "contrast", "id" => "g-refgroup_g-nonrefgroup"
			#    but: 
			#        "name" => "'non-ref' vs 'ref'" 
			} else {
				foreach my $arrOrgTime (sort keys %H_assayGroups1DifferenceTimeSeries) {
					foreach my $FV1 (sort keys %{$H_assayGroups1DifferenceTimeSeries{$arrOrgTime}}) {
					
						## Get FV2 assay_group
						my $fv1Group = $H_inAssayGroup{$arrOrgTime}{$FV1} ;
						foreach my $FV2 (sort keys %{$H_assayGroups1DifferenceTimeSeries{$arrOrgTime}{$FV1}}) {
			
							if ($arrOrgTime =~ /$array/) {

								## If pair never seen
								if (!exists $H_assayGroupDone{$arrOrgTime}{$FV1}{$FV2}) {
						
									## Store that pair so won't do it again (A,B or B,A)
									$H_assayGroupDone{$arrOrgTime}{$FV1}{$FV2} = 1 ;
									$H_assayGroupDone{$arrOrgTime}{$FV2}{$FV1} = 1 ;
						
									## Get FV2 assay_group
									my $fv2Group = $H_inAssayGroup{$arrOrgTime}{$FV2} ;

									$writer->startTag("contrast", "id" => "g${fv1Group}_g$fv2Group") ;
									if ($debug) { &log($logFileHandle, "[DEBUG] Printing <contrast> - $FV2 (g$fv2Group) vs. $FV1 (g${fv1Group}) - no ref.\n") ; }
																	
									## If 1 array, or RNA_seq
									# Remove the ;;;
									$FV1 =~ s/;;;/; /g ;
									$FV2 =~ s/;;;/; /g ;

									## Replace the \t\t separators by ;
									# ...so match the label in <assay_group id="gZ" label="xxx; yyy"> 
									$FV1 =~ s/\t\t/; /g ;
									$FV2 =~ s/\t\t/; /g ;
									if ($arrayNumber == 1) { $writer->dataElement("name" => "'$FV2' vs '$FV1'") ; }
		
									## If > 1 array,
									# user friendly name for the array
									else {
										if ($H_arrayIDs2arrayNames{$array} ne '') {
											$writer->dataElement("name" => "'$FV2' vs '$FV1' on '$H_arrayIDs2arrayNames{$array}'") ;
                    					} else { &log($logFileHandle, "[ERROR] $expAcc - No user-friendly name for array $array") ; die ; }
									}
									$writer->dataElement("reference_assay_group" => "g$fv1Group") ;
									$writer->dataElement("test_assay_group" => "g$fv2Group") ;
									$writer->endTag("contrast") ;
								}
							}
						}
					}
				}
			}
			if ($debug) { &log($logFileHandle, "[DEBUG] Closing <contrastS>\n") ;  }
			$writer->endTag("contrasts") ;
		}
		$writer->endTag("analytics") ;
	}
	$writer->endTag("configuration") ;
	$writer->end ;

	## Close file
	$XML->close ;

} else {
	&log($logFileHandle, "[WARNING] No config for $expAcc - not enough assay left after parsing\n") ;
}



## Subroutine
#############
## Print usage for the program
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
		"\t-debug: debug mode.\n") ;
}

## Log file
# Dual log to STDOUT and $logFileHandle. The latter will be reported by Conan in the interface; the former will 
# be reported in the LSF output file (and via an email from Conan) in the case of failure of this script.
sub log() {
	my $logFileHandle = shift;
	my $msg = shift;
	print STDOUT "$msg\n";
	printf($logFileHandle "$msg\n");
}

## Returns the number of values in a hash of hash
# $H{A}{A1}=un
# $H{A}{A2}=deux
# $H{B}{B1}=trois
# $H{B}{B2}=quatre
# -> returns 4
sub valuesNumberHashOfHash {
	my $hashOfHash = $_[0] ;
	my $nValues ;
	foreach my $k (keys %$hashOfHash) { $nValues += values %{$$hashOfHash{$k}} ; }
	return $nValues ;
}


## Compare 2 arrays
# (= identical values at the same position)
# and return the number of differences
sub compareLists {
	
	## Array1 is the one for which we return the results
	my $stringFV1 = $_[0] ;
    my $stringFV2 = $_[1] ;
	my $difference = 0 ;

	my @arrayFV1 = split("\t\t", $stringFV1) ;
    my @arrayFV2 = split("\t\t", $stringFV2) ;

	## Loop through array1, compare to array2, 
	# and count the differences
	foreach my $cpt (0..$#arrayFV1) {
		if ($arrayFV1[$cpt] ne $arrayFV2[$cpt]) {
			$difference++ ;
		}
	}

	## If array2 longer than array1, add the differences
	if (scalar @arrayFV2 > scalar @arrayFV1) { 
		$difference += ((scalar @arrayFV2) - (scalar @arrayFV2)) ; 
	}

	## Return number of differences
	return ($difference)	
}


## Read magetab file (SDRF / IDF files)
sub readMagetab {
    my $efvs2runAccessions = {} ;
	my %H_fv2fvType ;
	my $factorTypeString = "" ;
    my %H_technicalReplicateGroup ;

	my $i = -1 ; #counter for factor type (position) - -1, so array starts at 0
	my %H_factorList ;
	my %H_factorTypesOrder ; #record the factor type order
	my @A_factorValueInOrder ; #store factor values, in order
	my %H_assayFV ; #store the assay and the order of the FTypes

	if ($debug) { &log($logFileHandle, "\n[DEBUG] ===== READING MAGETAB (readMagetab module) =====\n") ; }

	## Create a Magetab4Atlas object. This reads the MAGETAB documents (via the
	# IDF filename pro vided) to get Atlas-relevant information for each assay in
	# the experiment.
	#
	# Catch the error from the magetab4atlas module
	# ... and only print the 1st line
	my $magetab4atlas ;
	eval { $magetab4atlas = Magetab4Atlas->new( "idf_filename" => $idf );  }; 
	if ($@)  { 
		my @A_moduleErrorMessage = split("\n", $@) ;
		&log($logFileHandle, "[ERROR] $expAcc -- Error in the Magetab files. $A_moduleErrorMessage[0]") ; die ;
	}

	## Test if we find any assays.
	if(!$magetab4atlas->has_assays) { &log($logFileHandle, "[DEBUG] No assay found!\n") ; } #die "No assay for this experiment!\n"; }
	 
	## Experiment type
	if ($debug) { &log($logFileHandle, "[DEBUG] Experiment type: ".$magetab4atlas->get_experiment_type."\n") ; }		
	my $expType = $magetab4atlas->get_experiment_type ; 
	
	my @A_magetabAssay = $magetab4atlas->get_assays ;
	if ($debug) { &log($logFileHandle, "[DEBUG] Assays: $A_magetabAssay[0]\n") ; }

	## Get the assay, or die
    if (!@{ $magetab4atlas->get_assays }) { &log($logFileHandle, "[ERROR] $expAcc - Cannot extract assay: no name or no factor values") ; die ; }

	foreach my $assay4atlas (@{ $magetab4atlas->get_assays }) {
		if ($debug) { &log($logFileHandle, "[DEBUG] Assays found !\n") ; }
		my @A_factorValueInOrder ; 

		## Get the organism
		my $organism = $assay4atlas->get_organism() ; 
        if ($debug) { &log($logFileHandle, "[DEBUG]\tOrganism:\n\t\t", $assay4atlas->get_organism, " ($organism) \n") ; }

		## Get the assay name
		#   Old SDRF: should be ENA_run (RNA-seq) or hybridisation name (microarray)
		#   Newer SDRF: assay name field
		#  For RNA-seq, also get the number of fastqURI - should be 1 (SE) or 2 (PE)
		my %H_runAccession ;
		if ($debug) { &log($logFileHandle, "[DEBUG] Assay: ", $assay4atlas->get_name, " (",$assay4atlas->get_name(),")\n") ; }	
		if ($expType eq "RNA-seq") {
			if($assay4atlas->has_fastq_uri_set) {
				foreach my $ENArun (keys %{ $assay4atlas->get_fastq_uri_set }) {
					foreach my $fastqUri (@{ $assay4atlas->get_fastq_uri_set->{ $ENArun } }) {
						$H_runAccession{$ENArun}++ ; #this is to record PE (2) or SE (1)
					}
					if ($debug) { &log($logFileHandle, "[DEBUG]\tENA run: $ENArun ($H_runAccession{$ENArun})\n") ; } 
				}
			}
		} else {
        	$H_runAccession{$assay4atlas->get_name()} = 1 ;
		}

		## Technical replicates, if any
		# They are linked to assay, but storing them associated to the run accession
		#  ...as they are going to be printed associated to the run accession
		#  (same tech. replicate group for all runs of an assay)
		if ($assay4atlas->has_technical_replicate_group) {
			my $technicalReplicateGroup = $assay4atlas->get_technical_replicate_group ;
			foreach my $run (keys %H_runAccession) {
				$technicalReplicateGroup =~ s/group/t/ ;
				$technicalReplicateGroup =~ s/\s//g ;
				$H_technicalReplicateGroup{$run} = $technicalReplicateGroup ;
				if ($debug) { &log($logFileHandle, "[DEBUG] $run storing technical replicate group $technicalReplicateGroup\n") ; }
			}
		}

		## Get the Factor type(s) and values for this assay
		if ($debug) { &log($logFileHandle, "[DEBUG]\tFactors:\n") ; }
		my $H_factors = $assay4atlas->get_factors;
		my $factorValueString = "" ; 
		my $factorTypeString = "" ;
		my $factorTime = "NOTIME" ;

		foreach my $factorType (keys %{ $H_factors }) {

			## If doesn't contain time as F.V.
			if (lc($factorType) !~ /time/) {

				## New f.v. string, '\t\t' separated
				$factorValueString .= $H_factors->{ $factorType }."\t\t" ;
            	$factorTypeString .= $factorType."\t\t" ;

				if ($debug) { &log($logFileHandle, "[DEBUG]\t\t$factorType: ", $H_factors->{ $factorType }, "\n") ; }

				## If 1st time: 
				# ... store the factor types, IN ORDER!
				# ... and the factor values in the same order
				## If not 1st time, 
				# ... store the factor values in the same order as the factor types 
				if (!exists $H_factorList{$factorType}) {

					## Check if ok to do it that way
					$i++ ;
					$H_factorList{$factorType} = $i ;
                	if ($debug) { &log($logFileHandle, "[DEBUG] New factor type '$factorType' - adding to the list position $i\n") ; }

					## ...and store
					$A_factorValueInOrder[$i] = $H_factors->{ $factorType } ;
					if ($debug) { &log($logFileHandle, "[DEBUG] Storing at position $i ".$H_factors->{ $factorType }."\n") ; }
				} else {

					## Position of that factor
					my $pos = $H_factorList{$factorType} ;
					
					## Store at that position
					$A_factorValueInOrder[$pos] = $H_factors->{ $factorType } ;
                	if ($debug) { &log($logFileHandle, "[DEBUG] (known factor) Storing at position $pos: ".$H_factors->{ $factorType }."\n") ; }
				}

			## If time factor, store it in a special variable
			} else {
				$factorTime = $H_factors->{ $factorType } ;
			}
		}

		## Store them associated to assay_name
		my $assayName = $assay4atlas->get_name() ;

		## Store them associated to array_design, 
		# Get it for a microarray assay
		my $arrayDesign = 0 ;
		if($expType =~ /array/) {
			if($assay4atlas->has_array_design) {
				$arrayDesign = $assay4atlas->get_array_design() ;
				if ($debug) { &log($logFileHandle, "\tArray design:\n\t\t", $assay4atlas->get_array_design, " ($arrayDesign) \n") ; }
			}
		}

		## String of array/organism/time
		my $arraySpeciesTime = "$arrayDesign;$organism;$factorTime" ;

		## Since storing a ref to the array doesn't work,
		# Make a string of the array, and store it
		my $ListFactorValueInOrder = join("\t\t", @A_factorValueInOrder);
		$H_assayFV{$arraySpeciesTime}{$assayName} = $ListFactorValueInOrder ;  

		$factorValueString =~ s/\t\t$// ; #remove trailing characters
		$factorTypeString =~ s/\t\t$// ; #remove trailing characters
		if ($debug) { &log($logFileHandle, "\t\tfactorValue string: $factorValueString\n") ; }
        if ($debug) { &log($logFileHandle, "\t\tfactorType string: $factorTypeString\n") ; }

		## For each factorValueString, store the factor types and their value
		# key is factorValueString because it needs to link to another hash which key is also factorValueString
		$H_fv2fvType{$arraySpeciesTime}{$ListFactorValueInOrder} = $H_factors ;

		## Store
		# Go through %H_runAccession to store all runAccession	
		# If PE/SE requirement, parse it here
		foreach my $runAcc (keys %H_runAccession) {
			if ( (!$pese) || ($pese eq "PE" && ($H_runAccession{$runAcc} == 2)) || ($pese eq "SE" && ($H_runAccession{$runAcc} == 1)) ) {

				if(exists($efvs2runAccessions->{ $arrayDesign }->{ $organism }->{ $factorTime }->{ $ListFactorValueInOrder })) {
					push @{ $efvs2runAccessions->{ $arrayDesign }->{ $organism }->{ $factorTime }->{ $ListFactorValueInOrder } }, $runAcc ;
				} else {
					$efvs2runAccessions->{ $arrayDesign }->{ $organism }->{ $factorTime }->{ $ListFactorValueInOrder } = [ $runAcc ] ;	
				}
			}	
		}
	}

	## Return values:
	# experiment type, factor value types (e.g. "genotype" - if multiple array: same F.V. types because it's the same SDRF file)
	# and the reference to a hash of mappings between factor values and run accession
	return($expType, $factorTypeString, $efvs2runAccessions, \%H_fv2fvType, \%H_technicalReplicateGroup, \%H_assayFV) ;
}
## End of readMagetab (TAG, leave this comment)
