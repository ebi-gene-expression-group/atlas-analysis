#!/bin/bash

# A wrapper script for the experiment summary
# "strict mode"

set -euo pipefail

if [ $# -lt 1 ]; then
        echo "Usage: $0 PATH_TO_EXPERIMENT "
        echo "e.g. $0 ${ATLAS_PROD}/analysis/baseline/rna-seq/experiments/E-MTAB-513"
        exit 1;
fi

expPath=$1
e=`basename ${expPath}`
output="$expPath/$e-atlasExperimentSummary.Rdata"

"$(dirname "${BASH_SOURCE[0]}")"/createAtlasExperimentSummary.R \
	--source "$expPath"\
	--accession "$e" \
	--output "$output"

test -s "$output" || (  >&2 echo "$0 Failed to create summary: $expPath" ; exit 1 )
