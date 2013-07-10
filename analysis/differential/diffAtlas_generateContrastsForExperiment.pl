#!/usr/bin/env perl

# POD documentation - main docs before the code
=pod

=head1 NAME

  Karyn Megy - 18-June-13
  kmegy@ebi.ac.uk
  Ereap_Create_ContrastFile.pl

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

  Ereap_Create_ContrastFile.pl -idf file.idf.txt -conf config_file.txt -out contrast_file.xml

  E.g. 
  /nfs/ma/home/kmegy/eREAP/ereap/scripts/Ereap_Create_ContrastFile.pl -sdrf /nfs/ma/home/kmegy/eREAP_test_files/E-MTAB-141.sdrf.txt -conf /nfs/ma/home/kmegy/eREAP_test_files/Contrast_config_file.txt

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
my ($idf, $conf, $xml, $help, $command_line) ; #arguments
my %H_Config ; #contains info. from the config file 


## Get arguments
################
GetOptions( 'help|Help|h|H' => \$help,
	    'idf=s' => \$idf,
            'out=s' => \$xml,
            'conf=s'=> \$conf
          ) ;

$command_line = join(' ',@ARGV); 

if (!$idf || !$conf) { print "[WARNING] Missing idf (-idf $idf) or configuration file (-conf $conf)!\n" ; $help  = 1 ; }

if ($help) { usage($command_line) ; exit ; }


## Extract information from the config file
## Only getting synonyms for 'reference' for now, but likely to expand.
###########################################
open (CONF, $conf) || die "Can't open the configuration file $conf!\n" ;
while (my $line=<CONF>) {
        chomp $line ;   
        #Getting synonyms for 'reference' - and store them
        if ($line =~ /^REFERENCE/) { 
                $line = lc($line);
                my ($KeY, $VaLue) = split ("\t", $line) ;
                my @A_VaLue = split (";", $VaLue) ;
                foreach my $i (0..$#A_VaLue) {
                        $A_VaLue[$i] =~ s/^ // ; #rm potential first/last spaces
                        $A_VaLue[$i] =~ s/ $// ;
                        $H_Config{REFERENCE}{$A_VaLue[$i]} = 1 ;
                }       
        }       
}
close CONF ;


##Print for a test
#foreach my $word (keys %{$H_Config{REFERENCE}}) { print "->$word.\n" ; } print "=====\n" ;
#exit ;


## Test Atlas eligibility for that experiment
## Call AEAtlas.pm
############################################## 


#------------------------------------------------------------------ START COPY / PASTE

my $reader_params;

# Checker will always perform basic validation checks
# when the MAGE-TAB files are parsed.
# Here we specify that it runs additional checks as required
my $check_sets = {
	'EBI::FGPT::CheckSet::AEAtlas' => 'ae_atlas_eligibility',
};

# Set up parser params 
$reader_params->{'check_sets'} = $check_sets;
$reader_params->{'skip_data_checks'} = '1';  #Skips checking for the presence of raw and processed data files
$reader_params->{idf} = $idf ;
$reader_params->{data_dir} = $dir ;

print STDOUT "[FLAG] reader_param\n" ;

# Call checker
my $checker = EBI::FGPT::Reader::MAGETAB->new($reader_params);
print STDOUT "[FLAG] checker\n" ;

$checker->parse();
print STDOUT "[FLAG] parse\n" ;

# Prints errors and warnings to SDTOUT (basic parsing + Atlas-specific checks)
$checker->print_checker_status();

print STDOUT "[FLAG] print_checker_status\n" ;


#Exit with 1 (failed) or 0 (passed)
if ($checker->has_errors) { print STDOUT "FAILED !\n" ; } ##exit 1 ; }
else { print STDOUT "PASSED!\n" ; } ##exit 0 ; }

#------------------------------------------------------------------ STOP STOP STOP


## If passed, 
## check for 3 replicates (FactorValues)
## then generate contrast file
##################################


## Check Factor Values
#foreach my $sdrf_row ($magetab->get_sdrfRows){
#
#	my @all_factor_values = $sdrf_row->get_factorValues;            
#
#        foreach my $factor_value($sdrf_row->get_factorValues) {                      
#		if ($factor_value->get_term) {   # anything except measurements
#			my $factor_name = $factor_value->get_term->get_category;
#                        my $factor_type = $factor_value->get_factor->get_factorType->get_value;
#                        my $factor_value_term = $factor_value->get_term->get_value;
#			print "NAME: $factor_name, TYPE: $factor_type, VALUE TERM: $factor_value_term\n" ;
#		}
#
#		#if ($factor_value->get_measurement) { # For measurements (anything which come with units), e.g. Age, Dose, Time, etc 
#                #	my $factor_type = $factor_value->get_factor->get_factorType->get_value;
#                #        my $factor_value_measurement = $factor_value->get_measurement->get_value;
#            	#}             
#	}
#}


## Subroutine
#############
sub usage {
    my ($command_line) = @_;
    
    print "Your command line was:\t".
	"$0 $command_line\n".
	"Compulsory parameters:\n".
	"\t-idf: IDF file name. Note that there there should be an SDRF file in the same folder - it'll be picked automatically up by the program (MAGETAB modules)\n".
	"\t-conf: configuration file name\n".
	"\t-xml: xml output file name\n" ;

}

