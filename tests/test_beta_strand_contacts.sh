#!/bin/bash

# Test script for beta strand contact parameters in RFdiffusion
# This script tests the new configuration parameters for beta strand contacts:
# - scaffoldguided.set_binder_beta_sheet_length
# - scaffoldguided.set_target_beta_sheet
# Using the umami.pdb structure with the beta sheet at A5-10 (1-indexed)

# Set paths
mFCGR4_PDB="/Users/matthewnemeth/Documents/1_Projects/RFdiffusion-beta/mFCGR4.pdb"
mFCGR4_SS="/Users/matthewnemeth/Documents/1_Projects/RFdiffusion-beta/mFCGR4_ss.pt"
mFCGR4_ADJ="/Users/matthewnemeth/Documents/1_Projects/RFdiffusion-beta/mFCGR4_adj.pt"
OUTPUT_DIR="/Users/matthewnemeth/Documents/1_Projects/RFdiffusion-beta/outputs"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Activate the RFdiffusion-mac conda environment
echo "Activating RFdiffusion-mac conda environment..."
eval "$(conda shell.bash hook)"
conda activate RFdiffusion-mac

echo "==== Running Beta Strand Contact Test ===="

# Run RFdiffusion with beta strand contact parameters
cd "$(dirname "$0")/.."

scripts/run_inference.py \
  inference.output_prefix="$OUTPUT_DIR/beta_strand_test" \
  inference.input_pdb="$mFCGR4_PDB" \
  inference.num_designs=1 \
  ppi.hotspot_res="[A5,A6,A7,A8,A9,A10]" \
  scaffoldguided.scaffoldguided=True \
  scaffoldguided.target_pdb=True \
  scaffoldguided.target_path="$mFCGR4_PDB" \
  scaffoldguided.set_binder_beta_sheet_length=6 \
  scaffoldguided.set_target_beta_sheet="[A5-10/0]" \
  scaffoldguided.flexible_beta_sheet=False \
  scaffoldguided.target_ss="$mFCGR4_SS" \
  scaffoldguided.target_adj="$mFCGR4_ADJ" \
  diffuser.partial_T=10 \
#   logging.save_ss_adj=True

# Check if the run was successful
if [ $? -eq 0 ]; then
  echo "\nBeta strand contact test completed successfully!"
  echo "Output files are in $OUTPUT_DIR/beta_strand_test"
  
  # List the generated files
  echo "\nGenerated files:"
  ls -l "$OUTPUT_DIR/beta_strand_test"*
else
  echo "\nBeta strand contact test failed."
fi

# Deactivate the conda environment
conda deactivate
