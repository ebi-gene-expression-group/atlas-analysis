#!/bin/bash
# A wrapper script for coexpressions
# Use TPM file, if does not exist try the FPKM file (we do not expect results to differ)

if [ $# -lt 1 ]; then
        echo "Usage: $0 PATH_TO_EXPERIMENT "
        echo "e.g. $0 ${ATLAS_PROD}/analysis/baseline/rna-seq/experiments/E-MTAB-513"
        exit 1;
fi
expPath=$1
e=`basename ${expPath}`

FPKMexpressionsFile="${expPath}/${e}-fpkms.tsv.undecorated.aggregated"
TPMexpressionsFile="${expPath}/${e}-tpms.tsv.undecorated.aggregated"

outputPath="${expPath}/${e}-coexpressions.tsv.gz"

if [ -s "${expPath}/$TPMexpressionsFile" ]; then
	$(dirname "${BASH_SOURCE[0]}")/run_coexpression_for_experiment.R "${expPath}/$TPMexpressionsFile" "$outputPath"
elif [ -s "${expPath}/$FPKMexpressionsFile" ]; then
	$(dirname "${BASH_SOURCE[0]}")/run_coexpression_for_experiment.R "${expPath}/$FPKMexpressionsFile" "$outputPath"
else
	>&2 echo "$0 data file not found in $1" ; exit 1
if
