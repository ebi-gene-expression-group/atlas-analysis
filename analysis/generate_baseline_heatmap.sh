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

expressionsFileTpm="${expPath}/${e}-tpms.tsv"
expressionsFileFpkm="${expPath}/${e}-fpkms.tsv"

configurationXml="${expPath}/${e}-configuration.xml"
outputPathTpm="${expPath}/${e}-heatmap-tpms.pdf"
outputPathFpkm="${expPath}/${e}-heatmap-fpkms.pdf"

if [ -s "$expressionsFileTpm" ]; then
	$(dirname "${BASH_SOURCE[0]}")/generateBaselineHeatmap.R --configuration "$configurationXml" \
		--input "$expressionsFileTpm" \
		--output "$outputPathTpm"
fi
if [ -s "$expressionsFileFpkm" ]; then
	$(dirname "${BASH_SOURCE[0]}")/generateBaselineHeatmap.R --configuration "$configurationXml" \
		--input "$expressionsFileFpkm" \
		--output "$outputPathFpkm"
fi

outputPath="${expPath}/${e}-heatmap.pdf"

if [-s "$outputPathTpm"]; then
	ln -s "$outputPathTpm" "$outputPath"
elif [-s "$outputPathFpkm"]; then
	ln -s "$outputPathFpkm" "$outputPath"
else
	>&2 echo "$0 No heatmap produced! $1" ; exit 1
fi
