#!/bin/bash    

source ${ATLAS_PROD}/sw/atlasprod/db/scripts/experiment_loading_routines.sh

# A wrapper script to allow for re-computing of differential expression for a microarraye experiment - so that this can be called from decorate_all_experiments.sh
expAcc=$1
expTargetDir=$2

if [ $# -lt 2 ]; then
  echo "Usage: $0 expAcc expPath"
  exit 1
fi 

# Calculate analytics
pushd ${expTargetDir}
diffAtlas_DE.pl -atlasXML=${expAcc}-configuration.xml
if [ $? -ne 0 ]; then
     echo "ERROR: Failed to calculate analytics for ${expAcc}" >&2
     exit 1
fi

# Generate experiment summary R object
createAtlasExperimentSummary.R $expAcc
if [ $? -ne 0 ]; then
   echo "ERROR: Failed to create an R experiment summary for ${expAcc}" >&2
   rm -rf *.Rdata
   exit 1
fi
popd

echo "Decorating .undecorated data files for experiment: $expAcc"
decorate_microarray_experiment.sh ${ATLAS_PROD}/analysis/differential/microarray/experiments/$expAcc
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to decorate ${ATLAS_PROD}/analysis/differential/microarray/experiments/$expAcc" >&2
    exit 1
fi

# Perform GSEA against Go/Interpro/Reactome
gxa_generate_gsea_for_differential_experiment.sh ${ATLAS_PROD}/analysis/differential/microarray/experiments/$expAcc
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to perform GSEA for ${expAcc}" >&2
    rm -rf ${ATLAS_PROD}/analysis/differential/microarray/experiments/$expAcc/*gsea*
    exit 1
fi

# Generate Genome browser tracks (bedGraph) for the experiment
gxa_generate_tracks_for_experiment.sh ${ATLAS_PROD}/analysis/differential/microarray/experiments/$expAcc
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to generate tracks for: ${ATLAS_PROD}/analysis/differential/microarray/experiments/$expAcc" >&2
    exit 1
fi

# Record if $expAcc is private - we need to know this 
# 1. to set appropriate unix permissions for the experiment directory under $ATLAS_EXPS, and 
# 2. to load the experiment with an appropriate private/public status into the staging instance
isPrivate=`is_exp_private $expAcc`
if [ -z "$isPrivate" ]; then 
  echo "ERROR: Failed to retrieve public/private status for $expAcc" >&2
  exit 1
fi

echo "Copying data to stage for differential microarray experiment: $expAcc" 
copy_experiment_data_to_stage $expAcc microarray differential all $isPrivate
if [ $? -ne 0 ]; then
    echo "ERROR: Command failed: 'copy_experiment_data_to_stage $expAcc microarray differential all $isPrivate'" >&2
    exit 1
fi 