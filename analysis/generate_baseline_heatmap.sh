#!/bin/bash
# A wrapper script for the baseline heatmap
# "strict mode"
set -euo pipefail

if [ $# -lt 1 ]; then
        echo "Usage: $0 PATH_TO_EXPERIMENT "
        echo "e.g. $0 ${ATLAS_PROD}/analysis/baseline/rna-seq/experiments/E-MTAB-513"
        exit 1;
fi
expPath=$1
e=`basename ${expPath}`

TPMexpressionsFile="${expPath}/${e}-tpms.tsv"
FPKMexpressionsFile="${expPath}/${e}-fpkms.tsv"

configurationXml="${expPath}/${e}-configuration.xml"
outputPath="${expPath}/${e}-heatmap.pdf"

if [ -s "$TPMexpressionsFile" ]; then
	$(dirname "${BASH_SOURCE[0]}")/generateBaselineHeatmap.R \
		--input "$TPMexpressionsFile" \
		--configuration "$configurationXml" \
		--output "$outputPath"
elif [ -s "$FPKMexpressionsFile" ]; then
	$(dirname "${BASH_SOURCE[0]}")/generateBaselineHeatmap.R \
		--input "$FPKMexpressionsFile" \
		--configuration "$configurationXml" \
		--output "$outputPath"
else
	>&2 echo "$0 data file not found in $1" ; exit 1
fi
