#!/bin/bash
# A wrapper script - so that this can be called from recalculate_all_experiments.sh
expTargetDir=$1

if [ $# -lt 1 ]; then
  echo "Usage: $0 expPath"
  exit 1
fi

pushd $expTargetDir

$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/calculate_percentile_ranks_for_exp.pl

popd
