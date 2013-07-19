#!/usr/bin/env perl

# POD documentation - main docs before the code
=pod

=head1 NAME

  Karyn Megy - 19-July-13
  kmegy@ebi.ac.uk
  diffAtlas_generateContrastsForExperiment_parallel.pl

=head1 SYNOPSIS
  
  Run the program diffAtlas_generateContrastsForExperiment_parallel.pl on a serie of experiments
  
=head1 DESCRIPTION

  Run the program diffAtlas_generateContrastsForExperiment_parallel.pl on a serie of experiments
  The experiment name is taken from the SDRF files, all in the same directory. 

  diffAtlas_generateContrastsForExperiment_parallel.pl:
  	generate a contrast file (XML) from a SDRF file (ArrayExpress format) and a config file. 

  Note that by a slight modifications of the program it would be easy to get the experiment name from another type of file. 

=head1 OPTIONS

  none

=head1 EXAMPLES

  diffAtlas_generateContrastsForExperiment_parallel.pl -dir <directory with expriment name> -out <directory for output files>

E.g.: 
  diffAtlas_generateContrastsForExperiment_parallel.pl -dir /ebi/microarray/home/mkeays/Atlas/differential/microarray/singleFactorContrasts/done/

=cut


use strict ;
use Getopt::Long qw(:config no_ignore_case);


## Declare $, @ and %
my ($help, $dir, $out) ; #arguments
my $commandLine ;


## Get option
GetOptions( 'help|Help|h|H' => \$help,
            'dir=s'  => \$dir,
          ) ;

$commandLine = join(' ',@ARGV); 

if (!$dir) { print "[WARNING] Missing directory (-dir)\n\n" ; $help  = 1 ; }

if ($help) { usage($commandLine) ; die ; }



## Get experiment SDRF files 
my $ls = `ls $dir/*sdrf*` ;
my @A_ls = split("\n", $ls) ;


## For each: extract the experiment name
## ... and run the program for contrast generation files 
foreach my $path (@A_ls) {
	my $experiment ;
	
	#Experiment name
	if ($path =~ /.+?\/(E-\w+?-\d+?).sdrf.txt$/) { $experiment = $1 ; }

	#Run program
	print "Experiment $experiment\n" ;
	my $CMD = "analysis/differential/diffAtlas_generateContrastsForExperiment.pl -exp $experiment  -conf analysis/differential/reference_assay_group_factor_values.txt" ;
	system($CMD) ;

	#Catch errors, if any
	if ($? != 0) { "[ERROR] Failed to execute $CMD. $!\n" ; }
	else { print  "Command $CMD successful\n" ; }
}



## Subroutine
#############
#Print usage for the program
sub usage {
	my ($command_line) = @_;

    	print "Your command line was:\t".
	    "$0 $command_line\n".
	    "Compulsory parameters:\n".
	    "\t-dir: directory containing the SDRF files for the experiments for which the XML contrast files will be generated.\n".
	    "Optional parameters:\n".
	    "\t-out: directory to write output files. Default is default output directory for diffAtlas_generateContrastsForExperiment_parallel.pl\n" ;
}



