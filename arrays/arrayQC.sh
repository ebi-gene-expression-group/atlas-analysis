#!/usr/bin/env bash

if [ $# -lt 1 ]; then
   echo "Usage: $0 expAcc"
   echo "e.g. $0 E-MTAB-1066 E-MTAB-1066.idf.txt AE_loadDirectoryPath "
   exit 1
fi

expTargetDir=$1
idf_filename=$2
ae_dir=$3
mirbase_dir=$4
scriptDir=$5

expAcc="$(basename $expTargetDir)"


rm -rf $expTargetDir/qc

#pushd $expTargetDir || exit 1 > /dev/null

perl $scriptDir/arrayQC.pl $expTargetDir $idf_filename $ae_dir $mirbase_dir
exitCode=$?
if [ $exitCode -eq 1 ]; then
    # The QC procedure succeeded but the experiment failed the QC
    #popd || exit 1 > /dev/null
    mv $expTargetDir ${ATLAS_PROD}/failedQC/microarray/
    echo "[QC] Quality control for ${expAcc} has failed - see ${ATLAS_PROD}/failedQC/microarray/${expAcc} for more info"
    exit 2
elif [ $exitCode -ne 0 ]; then
    #popd || exit 1 > /dev/null
    # The QC procedure itself failed (e.g. due to lack of memory) - we don't know if the experiment passes or fails the QC
    # Perl die() returns exit code=255
    echo "ERROR: Failed to perform QC for ${expTargetDir} - exit code: $exitCode" >&2
    exit 1
else
    # Experiment has passed QC check - move quality report dir into qc/
    find . -name "$expAcc*_QM" -type d | while read -r qcDir; do
        destination="$expTargetDir/qc/$(basename $qcDir)"
        mkdir -p $destination
        cp $qcDir/* $destination
        test -e "$qcDir/arrayQualityMetrics.js" \
            > "$destination/arrayQualityMetrics.js"
        rm -rf "$qcDir"
    done
    #popd || exit 1 > /dev/null
fi
