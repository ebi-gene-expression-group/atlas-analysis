#!/bin/bash

set -euo pipefail

#scriptDir=$2/bin  # $(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
#projectRoot=${scriptDir}/..

path_to_experiment_config=""
path_to_source=""
path_to_destination=""
usage(){
	[ ! "$1" ] || echo "${BASH_SOURCE[0]} : $1" >&2
	echo "${BASH_SOURCE[0]} Usage: -c <configuration.xml> -s <source tsv> -d <destination> -b <scriptdir>" >&2
	exit 2
}

while getopts ":c:s:d:b:" opt; do
  case $opt in
    c)
      path_to_experiment_config=$OPTARG;
      ;;
    s)
      path_to_source=$OPTARG;
      ;;
	  d)
	    path_to_destination=$OPTARG;
	  ;;
    b)
	    scriptDir=$OPTARG;
	  ;;
    ?)
      usage "Unknown option: $OPTARG"
    ;;
  esac
done

[ -f "$path_to_experiment_config" ] || usage "experiment configuration file not found: $path_to_experiment_config"
[ -f "$path_to_source" ] || usage "data file not found: $path_to_source"

mkdir -p "$(dirname $path_to_destination)"

tmp_files="$path_to_destination.$(date -u +"%Y-%m-%d_%H:%M").tmp"
trap 'rm -f ${tmp_files}.*' INT TERM EXIT

tmp_config_tsv=$tmp_files.config.tsv
#"$scriptDir/get_assay_groups.pl" "$path_to_experiment_config" > "$tmp_config_tsv"
"$scriptDir/atlas-bash-util/get_assay_groups.pl" "$path_to_experiment_config" > "$tmp_config_tsv"

columns_with_right_assays=$(join -1 2 -o 1.1 -2 1 <(head -n1 "$path_to_source" | tr $'\t' $'\n' | cat -n | sort -k2) <(cut -f 3 "$tmp_config_tsv" | sort -u) | tr $'\n' ',' | sed 's/,$//')

[ "$columns_with_right_assays" ] || usage "ERROR: $path_to_source headers have no assays from $path_to_experiment_config"

tmp_trimmed_source_tsv=$tmp_files.in.$(basename "$path_to_source" )

cut -f "1,$columns_with_right_assays" "$path_to_source" > "$tmp_trimmed_source_tsv"

tmp_out=$tmp_files.out.$(basename "$path_to_source" )

tmp_config_format_for_r=$tmp_files.config.tsv.manipulated
echo -e "AssayGroupID\tColumnHeading" > "$tmp_config_format_for_r"
cut -f 1,3 "$tmp_config_tsv" >> "$tmp_config_format_for_r"

Rscript "$scriptDir/atlas-analysis/norm/gxa_quantileNormalization.R" "$tmp_trimmed_source_tsv" "$tmp_config_format_for_r" "$tmp_out"

numRowsBeforeQuantileNormalization=$(wc -l < "$tmp_trimmed_source_tsv" )
numRowsAfterQuantileNormalization=$(wc -l < "$tmp_out" )
[ "$numRowsBeforeQuantileNormalization" -eq "$numRowsAfterQuantileNormalization" ] \
	|| usage "ERROR: attempted quantile normalization on file with $numRowsBeforeQuantileNormalization rows, got output with $numRowsAfterQuantileNormalization rows"

mv "$tmp_out" "$path_to_destination"

