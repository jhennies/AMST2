#!/bin/bash

# Convert a tif stack to ome-zarr (this is required for the pre-alignment workflow)
# For efficiency make sure the z-chunk size is 1!
snk_stack_to_ome_zarr \
  /path/to/tif_stack/ \
  /path/to/results/ \
  --batch_size 64 \
  --resolution 0.008 0.008 0.008 \
  --unit micrometer \
  --downsample_type Sample \
  --downsample_factors 2 2 2 \
  --chunk_size 1 512 512 \
  --cluster slurm \
  --cores 64 \
  --max_cores_per_task 16

# Run the alignment
snk_default_amst_pre_alignment \
  /scratch/hennies/tmp/amst2_test2/4t_raw.ome.zarr/ \
  /scratch/hennies/tmp/amst2_test2/ \
  --local_auto_mask \
  --template_roi 10.4 2.9 12.1 3.4 25.6 \
  --preview_downsample_level 3 \
  --batch_size 64 \
  --downsample_type Sample \
  --downsample_factors 2 2 2 2 2 \
  --chunk_size 1 512 512 \
  --cluster slurm \
  --cores 64 \
  --max_cores_per_task 16

# Optional: convert the result back to a tif stack
snk_ome_zarr_to_stack \
  /scratch/hennies/tmp/amst2_test2/pre-align.ome.zarr/ \
  /scratch/hennies/tmp/amst2_test2/ \
  pre-align-tiffs \
  --cores 64 \
  --batch_size 64 \
  --cluster slurm
