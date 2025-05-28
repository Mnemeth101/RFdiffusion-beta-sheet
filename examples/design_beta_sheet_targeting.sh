#!/bin/bash

# Set environment variables for your files and directories
OUTPUT_DIR="./outputs/beta_sheet_targeting"
mFCGR4_PDB="./inputs/your_target_structure.pdb"
mFCGR4_SS="./inputs/your_target_ss.pt" # You can generate this using the make_secstruc_adj.py script in the helper_scripts directory
mFCGR4_ADJ="./inputs/your_target_adj.pt" # You can generate this using the make_secstruc_adj.py script in the helper_scripts directory

# Make output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Run the beta-sheet targeting design
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
  logging.save_ss_adj=True
