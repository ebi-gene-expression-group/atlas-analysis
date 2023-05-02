#!/bin/bash
set -e

# Assign input parameters to variables
tissue=$1
accession=$2
sc_reference_C1=$3
sc_reference_C0=$4
sc_reference_phenData=$5
workflow_basedir=$6

# Set default status to success
export DECONV_STATUS="deconvolution successful"

mkdir -p scratch/${accession}/${accession}-${tissue}_scratch

# Run FARDEEP in the background
echo "start running FARDEEP"
fardeep_output=$(Rscript ${workflow_basedir}/atlas-analysis/deconvolution/FARDEEP_run.R \
  Tissue_splits/${accession}/${accession}-${tissue}-fpkms_scaled.rds $sc_reference_C1 \
  Output/${accession}/${accession}-${tissue}_res_FARDEEP.rds 2>&1) && { 
    echo "FARDEEP execution succeeded" 
} || { 
    echo "FARDEEP execution failed with error:" >&2 
    echo "$fardeep_output" >&2 
    export DECONV_STATUS="FARDEEP failed"
    exit 0
}

# Run DWLS in the background (that takes some time...)
echo "start running DWLS, that takes some time.."
dwl_output=$(Rscript ${workflow_basedir}/atlas-analysis/deconvolution/DWLS_run.R \
  Tissue_splits/${accession}/${accession}-${tissue}-fpkms_scaled.rds \
  $sc_reference_C0 $sc_reference_phenData 32 \
  Output/${accession}/${accession}-${tissue}_res_DWLS.rds \
  scratch/${accession}/${accession}-${tissue}_scratch 2>&1) && { 
    echo "DWLS execution succeeded" 
} || { 
    echo "DWLS execution failed with error:" >&2 
    echo "$dwl_output" >&2 
    export DECONV_STATUS="DWLS failed"
    exit 0
}

# Run EpiDISH in the background
echo "started running EpiDISH"
epidish_output=$(Rscript ${workflow_basedir}/atlas-analysis/deconvolution/EpiDISH_run.R \
  Tissue_splits/${accession}/${accession}-${tissue}-fpkms_scaled.rds \
  $sc_reference_C1 \
  Output/${accession}/${accession}-${tissue}_res_EpiDISH.rds 2>&1) && { 
    echo "EpiDISH execution succeeded" 
} || { 
    echo "EpiDISH execution failed with error:" >&2 
    echo "$epidish_output" >&2 
    export DECONV_STATUS="EpiDISH failed"
    exit 0
}

# Wait for all background processes to finish
wait

echo "Deconvolution status: $DECONV_STATUS"
