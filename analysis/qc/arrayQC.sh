#!/usr/bin/env bash

# Source script from the same (prod or test) Atlas environment as this script
scriptDir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source ${scriptDir}/../../bash_util/generic_routines.sh
atlasEnv=`atlas_env`

source ${ATLAS_PROD}/sw/atlasinstall_${atlasEnv}/atlasprod/irap/gxa_preirap.conf

if [ $# -lt 1 ]; then
   echo "Usage: $0 expAcc"
   echo "e.g. $0 E-MTAB-1066"
   exit 1
fi

expAcc=$1
pushd ${ATLAS_PROD}/analysis/differential/microarray/experiments/$expAcc
rm -rf qc

${ATLAS_PROD}/sw/atlasinstall_${atlasEnv}/atlasprod/analysis/qc/arrayQC.pl $expAcc
exitCode=$?
if [ $exitCode -eq 1 ]; then
    # The QC procedure succeeded but the experiment failed the QC
    popd
    mv ${ATLAS_PROD}/analysis/differential/microarray/experiments/$expAcc ${ATLAS_PROD}/failedQC/microarray/
    echo "[QC] Quality control for ${expAcc} has failed - see http://www.ebi.ac.uk/~rpetry/atlas3/failedQC/microarray/${expAcc} for more info"
    exit 2
elif [ $exitCode -ne 0 ]; then
    popd
    # The QC procedure itself failed (e.g. due to lack of memory) - we don't know if the experiment passes or fails the QC
    # Perl die() returns exit code=255
    echo "ERROR: Failed to perform QC for ${expAcc} - exit code: $exitCode" >&2
    exit 1       
else
    # Experiment has passed QC check - move quality report dir into qc/
    mkdir -p qc
    mv *_QM qc
    popd
fi