# Pre-alignment and AMST

This workflow contains of the following basic steps:

 1. Conversion to ome.zarr
 2. Template matching alignment or long-range slice-by-slice alignment
 3. Stack alignment with Simple Elastix
 4. Pre-alignment by combining Template matching (step 2) and Elastix stack alignment (step 3)
 5. Alignment to median smoothed template (AMST)

## Pre-requesites

Make sure the package is properly installed (see [README.md](../README.md)) 
and the AMST2 conda environment is active.

If this is the case, the following commands should be available and yield their respective help output:

```
amst-nsbs-pre-align -h
amst-cleanup-nsbs-pre-align -h
amst-run -h
amst-cleanup-run -h
```

## Pre-alignment (Steps 1 to 4)

Set up an output directory and cd into it

```
mkdir pre-align-out-dir
cd pre-align-out-dir
```

Copy and modify the [pre-alignment parameter file](../parameter_files/params-nsbs-pre-align-local.yaml)

```
cp /path/to/this/repo/parameter_files/params-nsbs-pre-align-local.yaml ./my-pre-align-params.yaml
vi ./my-pre-align-params.yaml
```

For the example dataset, only modify the input and output locations.

Run the pre-alignment:

```
amst-nsbs-pre-align ./my-pre-align-params.yaml
```

Once it's done, you can check the success by one of the following ways:

 1. Load the full-sized tif stack result stored in ```pre_align_to_tif_stack/nsbs-pre-align``` into Fiji (preferentially as virtual stack)
 2. Load the full-sized ome.zarr result ```hela-1-nsbs/apply_pre_align/nsbs-pre-align.ome.zarr``` with an ome.zarr viewer (e.g. BigDataViewer in Fiji)
 3. Load one of the previews of the intermediate steps ```hela-1-nsbs/sbs_alignment/sbs.meta/elastix_preview.h5``` (e.g. with Fiji -> File -> Import -> Hdf5...)

## AMST (step 5)

Set up an output directory and cd into it

```
mkdir amst-out-dir
cd amst-out-dir
```

Copy and modify the [AMST parameter file](../parameter_files/params-amst-bspline-local.yaml)

```
cp /path/to/this/repo/parameter_files/params-amst-bspline-local.yaml ./my-amst-params.yaml
vi ./my-amst-params.yaml
```

Point the inputs to the previously generated pre-alignment outputs which can be found here:

```
  input_dirpath: pre-align-out-dir/stack_to_ome_zarr/input-raw.ome.zarr
  pre_align_dirpath: pre-align-out-dir/apply_pre_align/nsbs-pre-align.ome.zarr
  pre_align_transforms: pre-align-out-dir/nsbs-pre-align.json
```

Run AMST:

```
amst-run ./my-amst-params.yaml
```

Now you can check the result:
 1. Full size tif-stack in ```hela-1-amst/amst_to_tif_stack/amst/```
 2. ome.zarr in ```hela-1-amst/apply_amst/amst.ome.zarr```

## Clean-up of intermediate results

To simplify the output folder structure as well as remove intermediate steps, run the following commands for 
pre-alignment and AMST individually:

```
amst-cleanup-nsbs-pre-align pre-align-out-dir/my-pre-align-params.yaml
amst-cleanup-run amst-out-dir/my-amst-params.yaml
```

Note that the directory structure changes such that re-running the AMST step would require adjusting the input 
file/directory locations.
