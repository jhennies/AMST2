#!/bin/bash

# Convert a tif stack to ome-zarr (this is required for the pre-alignment workflow)
# For efficiency make sure the z-chunk size is 1!
snk_stack_to_ome_zarr \
  /path/to/tif_stack/ \
  /path/to/results/ \
  --cores 64 \
  --batch_size 64 \
  --resolution 0.008 0.008 0.008 \
  --unit micrometer \
  --downsample_type Sample \
  --downsample_factors 2 2 2 \
  --chunk_size 1 512 512
  --cluster slurm \

# Run the alignment
snk_default_amst_pre_alignment \
  /scratch/hennies/tmp/amst2_test/subset_6318.ome.zarr/ \
  /scratch/hennies/tmp/amst2_test/ \
  --batch_size 64 \
  --cores 64 \
  --max_cores_per_task 16 \
  --local_auto_mask \
  --preview_downsample_level 3 \
  --template_roi 10.4 4.6 12.0 5.2 0.3 \
  --chunk_size 1 512 512 \
  --cluster slurm
