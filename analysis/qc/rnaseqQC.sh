#!/usr/bin/env bash

# Source script from the same (prod or test) Atlas environment as this script
scriptDir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source ${scriptDir}/../../bash_util/generic_routines.sh
atlasEnv=`atlas_env`

if [ $# -lt 1 ]; then
   echo "Usage: $0 expAcc"
   echo "e.g. $0 E-MTAB-1066"
   exit 1
fi

expAcc=$1
expTargetDir=$2
pushd ${expTargetDir}
rm -rf qc

${ATLAS_PROD}/sw/atlasinstall_${atlasEnv}/atlasprod/analysis/qc/rnaseqQC.pl $expAcc $expTargetDir > ${expAcc}.qc.log
exitCode=$?
if [ $exitCode -eq 1 ]; then
    rm -rf ${expAcc}.qc.log
	# The QC procedure succeeded but the experiment failed QC.
	popd
	mv ${expTargetDir} ${ATLAS_PROD}/failedQC/rna-seq
	echo "[QC] Quality control for ${expAcc} has failed - see http://www.ebi.ac.uk/~rpetry/atlas3/failedQC/rna-seq/${expAcc} for more info"
	exit 2
elif [ $exitCode -ne 0 ]; then
    rm -rf ${expAcc}.qc.log
	popd
	# The QC script itself failed (e.g. ran out of memory) -- we don't know
	# whether the experiment has passed or failed QC.
	echo "ERROR: Failed to perform QC for ${expAcc} -- exit code: $exitCode" >&2
	exit 1
else
	# Experiment has passed QC check -- 

    # See if there were any warnings about runs with less than 70% reads mapped.
    lowMQrunsFound=`grep LOW_MQ ${expAcc}.qc.log`

    if [ ! -z "$lowMQrunsFound" ]; then
        
# Use newline characters to separate elements in grep results.
IFS="
"
        for line in $lowMQrunsFound; do
            qcLine=${line/WARN  - /}
            echo $qcLine
        done
    fi

    # See if there are any warnings about QC fails.
    failedQCruns=`grep QC_FAIL ${expAcc}.qc.log`

    if [ ! -z "$failedQCruns" ]; then
IFS="
"
        for line in $failedQCruns; do
            failedLine=${line/WARN  - /}
            failedLine=${line/ERROR - /}
            echo $failedLine
        done
    fi

    # Print the info about percentage of runs that passed QC.
    pctPassed=`grep PCT_PASSED ${expAcc}.qc.log`

    if [ ! -z "$pctPassed" ]; then
IFS="
"
        for line in $pctPassed; do
            pctLine=${line/INFO  - /}
            echo $pctLine
        done
    fi
    
    rm -rf ${expAcc}.qc.log

    # Move the API results file to the qc/ dir.
	mkdir -p qc
	mv ${expAcc}-irap-single-lib-report.tsv qc
	popd
fi

