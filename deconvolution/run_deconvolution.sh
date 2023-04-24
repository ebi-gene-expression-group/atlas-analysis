#!/bin/bash

# exit when any command fails
set -e 
# Assign input parameters to variables
tissue=$1
accession=$2
sc_reference_C1=$3
sc_reference_C0=$4
sc_reference_phenData=$5
workflow_basedir=$6

# Run FARDEEP
Rscript ${workflow_basedir}/atlas-analysis/deconvolution/FARDEEP_run.R \
  Tissue_splits/${accession}/${accession}-${tissue}-fpkms_scaled.rds $sc_reference_C1 \
  Output/${accession}-${tissue}_res_FARDEEP.rds &
mkdir -p scratch/${accession}/${accession}-${tissue}_scratch &
  
# Run DWL (that takes some time...)
Rscript ${workflow_basedir}/atlas-analysis/deconvolution/DWLS_run.R \
  Tissue_splits/${accession}/${accession}-${tissue}-fpkms_scaled.rds \
  $sc_reference_C0 $sc_reference_phenData 32 \
  Output/${accession}-${tissue}_res_DWLS.rds \
  scratch/${accession}/${accession}-${tissue}_scratch &
  
# Run EpiDISH
Rscript ${workflow_basedir}/atlas-analysis/deconvolution/EpiDISH_run.R \
  Tissue_splits/${accession}/${accession}-${tissue}-fpkms_scaled.rds \
  $sc_reference_C1 \
  Output/${accession}-${tissue}_res_EpiDISH.rds

wait
