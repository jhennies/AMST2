#!/bin/bash

# The goal of this workflow is to introduce a long-range alignment (d=16 slices) in order to minimize the overall drift
# of the alignment.
# For all functions add parameters as required!

# Convert a tif stack to ome-zarr (this is required for the pre-alignment workflow)
# For efficiency make sure the z-chunk size is 1!
snk_stack_to_ome_zarr \
  /path/to/tif_stack/input_stack/ \
  /path/to/results/ \
  --resolution 0.008 0.008 0.008 \
  --unit micrometer \
  --downsample_type Sample \
  --downsample_factors 2 2 2 \
  --chunk_size 1 512 512

# Run alignment for local correspondences
# Using auto-pad here to make sure the data does not drift out of the canvas
snk_elastix_stack_alignment \
  /path/to/results/input_stack.ome.zarr \
  /path/to/results/ \
  --auto_mask \
  -out-oz-fn elastix-step1.ome.zarr \
  --apply_final \
  --downsample_type Sample \
  --auto_pad

# Run alignment for long-range correspondences
# This usually has much less bias to drift into a certain direction
# Switch of fixing of big jumps as we expect everything to be roughly in place at this stage
# No auto-pad as this does not work for z_step > 1, also we don't need it here
# Also, we don't need to apply the final alignment
snk_elastix_stack_alignment \
  /path/to/results/elastix-step1.ome.zarr \
  /path/to/results/ \
  --auto_mask \
  -out-oz-fn elastix-step16.ome.zarr \
  --z_step 16 \
  --downsample_type Sample \
  --no_fixing_of_big_jumps

# Join the transformations of the local and long-range alignment steps
# For this function, we don't need a snakemake workflow management, it's not computationally expensive,
#   so we just use the basic squirrel functionality
# Note that we have to keep the meta-data of the step1 alignment as  this contains the information of the bounding boxes
#   within each of the slices. This is needed for the auto-pad
linalg_dot_product_on_affines \
  /path/to/results/elastix-step1.meta/elastix.json \
  /path/to/results/elastix-step16.meta/elastix.json \
  /path/to/results/combined_alignment.json \
  --keep_meta 0
apply_auto_pad \
  /path/to/results/combined_alignment.json \
  /path/to/results/combined_alignment_autopad.json

# Apply the alignment
snk_apply_transformation \
  /path/to/results/input_stack.ome.zarr \
  /path/to/results/combined_alignment_autopad.json \
  /path/to/results/

# Optional: convert the result back to a tif stack
snk_ome_zarr_to_stack \
  /path/to/results/combined_alignment_autopad.ome.zarr \
  /path/to/results/ \
  combined_alignment_autopad
