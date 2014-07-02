#!/bin/bash    

# A wrapper script to allow for re-computing of differential expression for a microarraye experiment - so that this can be called from decorate_all_experiments.sh
$expAcc=$1
$expTargetDir=$2

if [ $# -lt 2 ]; then
  echo "Usage: $0 expAcc expPath"
  exit 1
fi 

pushd ${expTargetDir}
diffAtlas_DE.pl -atlasXML=${expAcc}-configuration.xml
if [ $? -ne 0 ]; then
     echo "ERROR: Failed to calculate analytics for ${expAcc}" >&2
     rm -rf *.tsv.undecorated
     rm -rf *.png
     exit 1
fi
popd