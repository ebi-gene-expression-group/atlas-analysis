#!/usr/bin/env perl

# POD documentation - main docs before the code
=pod

=head1 NAME

  Karyn Megy - 18-June-13
  kmegy@ebi.ac.uk
  diffAtlas_generateContrastsForExperiment.pl

=head1 SYNOPSIS

  Generate a contrast file (XML) from SDRF & IDF files (ArrayExpress format)

=head1 DESCRIPTION

  Generate a contrast file (XML) from a IDF & SDRF files (ArrayExpress format)
  The IDF file is the argument but both files, in the same folder, are necessary                                  
  The MAGETAB module will infer the SDRF file name from the IDF file.  

  Step-1: test Atlas eligibility 
  Step-2: if passed, test replicates, then generate contrast file

=head1 OPTIONS

  none

=head1 EXAMPLES

  gxa_ereap_create_contrastfile.pl -idf file.idf.txt -dir path/to/idf-sdrf-files -conf config_file.txt -out contrast_file.xml

  E.g. 
  diffAtlas_generateContrastsForExperiment.pl -idf E-MTAB-533.idf.txt -dir /net/isilon5/ma/home/arrayexpress/ae2_production/data/EXPERIMENT/MTAB/E-MTAB-533/ -conf atlasprod/analysis/differential/reference_assay_group_factor_values.txt

=cut

##Use specific PERL path
# --> future development: should we make it an option (centOS or redhat)?
#
#Running on banana (CentOS)
#use lib '/ebi/microarray/home/fgpt/sw/lib/perl/CentOS_prod/lib64/perl5/site_perl','/ebi/microarray/home/fgpt/sw/lib/perl/CentOS_prod/lib';
#use lib '/ebi/microarray/home/fgpt/sw/lib/perl/FGPT_CentOS_prod/lib';  # This takes precedence

#Running on ebi-001 (RedHat OS)
use lib '/ebi/microarray/home/fgpt/sw/lib/perl/RH_prod/lib64/perl5','/ebi/microarray/home/fgpt/sw/lib/perl/RH_prod/lib';
use lib '/ebi/microarray/home/fgpt/sw/lib/perl/RH_prod/lib';
use lib '/ebi/microarray/home/fgpt/sw/lib/perl/FGPT_RH_prod/lib/' ;
use lib '/nfs/ma/home/kmegy/' ;

use strict ;
use Getopt::Long qw(:config no_ignore_case);
use EBI::FGPT::Reader::MAGETAB;
use Bio::MAGETAB::Util::Reader;

use File::Spec; #required (?) for AEAtlas.pm
use Data::Dumper; #required (?) for AEAtlas.pm
use EBI::FGPT::Reader::MAGETAB; #required (?) for AEAtlas.pm
use EBI::FGPT::Common qw(date_now); #required (?) for AEAtlas.pm


# The following oracle-specific environment variables need to be set here at runtime
# or else all oracle DB connections (e.g. to Conan DB) needed in checks will fail.
BEGIN {
	$ENV{ORACLE_HOME}       = "/sw/arch/dbtools/oracle/product/9.2.0";	
	$ENV{ORA_NLS33}         = "/sw/arch/dbtools/oracle/product/9.2.0/ocommon/nls/admin/data";
}


## Initialise global $, @ and %
my ($idf, $dir, $conf, $xml, $help, $command_line) ; #arguments
my %H_Config ; #contains info. from the config file 
my %H_Concatenate_FactorValues ; #contains factor values to be parsed
my $flag ; #to parse config file
my $reference ; #reference factor value to calculate D.E. against
my $exit_code1 ; #exist code for Atlas checks
my $exit_code2 ; #exit code for FactorValue checks 


## Get arguments
################
GetOptions( 'help|Help|h|H' => \$help,
	    'idf=s' => \$idf,
            'out=s' => \$xml,	
            'dir=s' => \$dir,   
            'conf=s'=> \$conf
          ) ;

$command_line = join(' ',@ARGV); 

if (!$idf || !$conf || !$dir) { print "[WARNING] Missing directory (-dir $dir), idf (-idf $idf) or configuration file (-conf $conf)!\n" ; $help  = 1 ; }

if ($help) { usage($command_line) ; exit ; }

$idf = $dir."/".$idf ;

## Extract information from the config file
## Only getting synonyms for 'reference' for now, but likely to expand.
###########################################
open (CONF, $conf) || die "Can't open the configuration file $conf!\n" ;
while (my $line=<CONF>) {
        chomp $line ; 

	if ($line !~ /^#/) {
		if ($line =~ /REFERENCE/) { $flag = "REFERENCE" ; }
		elsif ($line =~ /FACTOR_VALUE_KILL/) { $flag = "FACTOR_VALUE_KILL" ; }
		elsif ($line =~ /\[\//) { $flag = "" ; }
        	else { if (($flag ne "") && ($line ne "")) { $H_Config{$flag}{$line} = 1 ; } }
	}
}	
close CONF ;


##Print for a test
#foreach my $CaT (keys %H_Config) {
#	print ">>$CaT<<\n" ; 
#	foreach my $VaLuE (keys %{$H_Config{$CaT}}) {
#		print ".$VaLuE.\n" ; 
#	}	
#} print "=====\n" ;
#exit ;


## Test Atlas eligibility for that experiment
## Call /ebi/microarray/home/fgpt/sw/lib/perl/FGPT_RH_prod/lib//EBI/FGPT/CheckSet/AEAtlas.pm
## Output file (ae_atlas_eligibility*) in the same directory as the IDF/SDRF files 
############################################## 
my $reader_params ;

# Checker will always perform basic validation checks
# when the MAGE-TAB files are parsed.
# Here we specify that it runs additional checks as required
my $check_sets = {
       'eREAP::ereap::scripts::AEAtlas' => 'ae_atlas_eligibility',
};

# Set up parser params 
$reader_params->{'check_sets'} = $check_sets;
$reader_params->{'skip_data_checks'} = '1';  #Skips checking for the presence of raw and processed data files
$reader_params->{idf} = $idf ;
$reader_params->{data_dir} = $dir ;

# Call checker

print STDOUT "[FLAG] Reader param: $check_sets ; $idf ; $dir\n" ; 

my $magetab = EBI::FGPT::Reader::MAGETAB->new($reader_params);

print STDOUT "[INFO] Checking Atlas eligibility\n" ;
$magetab->parse();

# Prints errors and warnings to STDOUT (basic parsing + Atlas-specific checks)
$magetab->print_checker_status();

#Return exit code: 1 (failed) or 0 (passed)
if ($magetab->has_errors) { print STDOUT "[INFO][CHECK] Atlas checks FAILED!\n" ; $exit_code1 = 1 ; } 
else { print STDOUT "[INFO][CHECK] Atlas checks PASSED!\n" ; }


## If Atlas checked passed (exit code 0), 
## Do more test on Factor Values (>= 3 replicates, reference etc.)
##################################
print STDOUT "[INFO] Gathering factor values\n" ;

print STDOUT "[FLAG] MAGETAB fetching!\n" ;
my $sdrf = $magetab->get_magetab ;
print STDOUT "[FLAG] MAGETAB fetched! $sdrf\n" ;
my $error_msge ;

$exit_code1 = "" ; #for testing
if (!$exit_code1) {

	## Check Factor Values
	foreach my $sdrf_row ($sdrf->get_sdrfRows){
		my @all_factor_values = $sdrf_row->get_factorValues;            
		my $factor_value_string ;

        	foreach my $factor_value($sdrf_row->get_factorValues) {                      

			if ($factor_value->get_term) {   # anything except measurements
				my $factor_name = $factor_value->get_term->get_category ;
                        	my $factor_type = $factor_value->get_factor->get_factorType->get_value ;
                        	my $factor_value_term = $factor_value->get_term->get_value ;
				print "$factor_name $factor_type $factor_value_term\n" ;	
	
				#Term belongs to the FactorValues kill list?
				foreach my $ref_value (keys %{$H_Config{FACTOR_VALUE_KILL}}) {
					if ($factor_type =~ /$ref_value/) { 
						$exit_code2 = 1 ; 
						$error_msge = "Experiment contains forbbiden term '$ref_value'\n" ;
						print "[INFO][CHECK] Factor Value check FAILED! $error_msge\n" ;
						exit 1 ;	
					}	
					
				}	
				
				#Concatenate the FactorValues
                                $factor_value_string = $factor_value_term ;
			} 

			if ($factor_value->get_measurement) { # For measurements (anything which come with units), e.g. Age, Dose, Time, etc 
                		my $factor_type = $factor_value->get_factor->get_factorType->get_value ;
                	      	my $factor_value_measurement = $factor_value->get_measurement->get_value ;
				my $factor_value_measurement_unit = $factor_value->get_measurement->get_unit->get_value ;
				#my $factor_value_measurement_unit_cat = $factor_value->get_measurement->get_unit->get_category ; #e.g. TimeUnit
				#print "$factor_type: $factor_value_measurement $factor_value_measurement_unit\n" ;
 
				#Concatenate the FactorValues
				$factor_value_string .= "-".$factor_value_measurement."-".$factor_value_measurement_unit ;
           		}      

			NODE: foreach my $node ($sdrf_row->get_nodes) {
				print "$node\n" ;
                        	if ( $node->isa('Bio::MAGETAB::Comment') ) {
                        		my $comment = $node;
                       			print "[COMMENT] ".$comment->get_name.": ".$comment->get_value."\n" ;
 				} 
			}

#			foreach my $comment ($sdrf->get_comments) {
#                        	if ($comment->get_name eq "ENA_RUN") {
#                                	print "[COMMENT] ".$comment->get_name.": ".$comment->get_value."\n" ;
#                        	}
#			}	
		}
		## Store the concatenated FactorValues
		$H_Concatenate_FactorValues{$factor_value_string}++ ;
                #print "Factor Value String: $factor_value_string\n" ;

		
		##Get the ENA run identifier
#		foreach my $comment ($sdrf->get_comments) {
#			if ($comment->get_name eq "ENA_RUN") {
#				print "[COMMENT] ".$comment->get_name.": ".$comment->get_value."\n" ;
#			}
#		}

	}
}	

## Parse the factor values: >= 3 replicates? Reference values?  
print "[INFO] Checking Factor Values suitability for differential expression analysis\n" ;
foreach my $FV_string (keys %H_Concatenate_FactorValues) {
	
	print "$FV_string\t$H_Concatenate_FactorValues{$FV_string}\n" ;

	#Test replicates, delete if < 3
	if ($H_Concatenate_FactorValues{$FV_string} < 3) { delete $H_Concatenate_FactorValues{$FV_string} ; }

	#Test reference	
	#Based on previous test, should have >=3 replicates!
	foreach my $ref_value (keys %{$H_Config{REFERENCE}}) {
		if ($FV_string =~ /$ref_value/) { 
			if (!defined $reference) { $reference = $FV_string ; }
			else { print STDOUT "MULTIPLE REFERENCE! $FV_string & $reference\n" ; } 
		}
	}	
}


#Reference Factor Value ? 
if (!defined $reference) { print STDOUT "[INFO][CHECK] Factor Value check FAILED! No reference ($reference)\n" ; $exit_code2 = 1 ; }

#Any Factor value left (>= 3 replicates)?
#Need at least 2 of them!
if (scalar %H_Concatenate_FactorValues < 2) { print STDOUT "[INFO][CHECK] Factor Value check FAILED! Less than 2 values with at least 3 replicates!\n" ; $exit_code2 = 1 ; }

#If no error reported, then FactorValue test passed!
if (!defined $exit_code2) { print "[INFO][CHECK] FactorValue test passed!\n" ;}

##If test passed, then format in XML
#print "<configuration>\n" ;
#print "\t<analytics>\n" ;

#print "\t\t<assay_groups>\n" ;
#....
#....
#print \t\t<\assay_groups>\n" ;


#print "\t\t<contrasts>\n" ;
#....
#....
#print \t\t<\contrasts>\n" ;


#print "\t<\analytics\n>"
#print "<\configuration>\n" ;


## Subroutine
#############
sub usage {
    my ($command_line) = @_;
    
    print "Your command line was:\t".
	"$0 $command_line\n".
	"Compulsory parameters:\n".
	"\t-idf: IDF file name [NO path!]\n".
        "\t-dir: directory in which the IDF and SDRF files are. The SDRF file will be picked up automatically by the program (MAGETAB modules)\n".
	"\t-conf: configuration file name\n".
	"\t-xml: xml output file name\n" ;

}

