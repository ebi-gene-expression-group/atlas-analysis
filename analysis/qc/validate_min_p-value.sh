#!/bin/bash
# Confirm that a tsv file has at least one p-value < 0.5
# If not it's not good enough for Atlas - exit with code 2 to signify this
# "strict mode"
set -euo pipefail

if [ $# -lt 1 ]; then
        echo "Usage: $0 PATH_TO_TSV " >&2
        echo "e.g. $0 $ATLAS_PROD/analysis/differential/microarray/experiments/E-ATMX-20/E-ATMX-20_A-AFFY-2-analytics.tsv" >&2
        exit 1
fi
path=$1

pValueColumns=$(head -n1 $path | tr $'\t' $'\n' | grep -n p-value | cut -f 1 -d : | tr $'\n' , | sed 's/,$//')
if [ ! "$pValueColumns" ]; then
	echo "ERROR: could not parse header (?)" >&2
	exit 1
fi

#pick columns with pvalue entries
#turn the columns into one long list
#skip NA values
#skip empty fields if any (shouldn't happen)
errorCode=2
cut -f $pValueColumns $path \
	| tail -n +1 \
	| tr $'\t' $'\n' \
	| grep -v NA \
	| grep -v -e '^[[:space:]]*$' \
	| while read -r field; do
	if [[ $(bc <<< "$field < 0.5" ) -eq 1 ]]; then
		#we're done
		errorCode=0
		break
	fi
done

exit $errorCode
