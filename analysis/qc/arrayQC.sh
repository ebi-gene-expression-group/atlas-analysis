#!/usr/bin/env bash

# Source script from the same (prod or test) Atlas environment as this script
scriptDir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
projectRoot=${scriptDir}/../..

if [ $# -lt 1 ]; then
   echo "Usage: $0 expAcc"
   echo "e.g. $0 E-MTAB-1066"
   exit 1
fi

expTargetDir=$1
expAcc="$(basename $expTargetDir)"

fixArrayQualityMetricsFile(){
	generatedFile=$1
	correctFile=$2
	matchPhrase='var highlightInitial\|var arrayMetadata\|var svgObjectNames'

	echo "/* patch part 1 - lines from generated file */"
	grep "$matchPhrase" < "$generatedFile"
    echo "/* patch part 2 - lines from patch file */"
	grep -v "$matchPhrase" < "$correctFile"
}

rm -rf $expTargetDir/qc

pushd $expTargetDir > /dev/null
$projectRoot/analysis/qc/arrayQC.pl $expAcc
exitCode=$?
if [ $exitCode -eq 1 ]; then
    # The QC procedure succeeded but the experiment failed the QC
    popd > /dev/null
    mv $expTargetDir ${ATLAS_PROD}/failedQC/microarray/
    echo "[QC] Quality control for ${expAcc} has failed - see http://www.ebi.ac.uk/~rpetry/atlas3/failedQC/microarray/${expAcc} for more info"
    exit 2
elif [ $exitCode -ne 0 ]; then
    popd > /dev/null
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
            && fixArrayQualityMetricsFile "$qcDir/arrayQualityMetrics.js" "$scriptDir/arrayQualityMetrics.js" \
            > "$destination/arrayQualityMetrics.js"
        rm -rf "$qcDir"
    done
    popd > /dev/null
fi
