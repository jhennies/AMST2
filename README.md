# AMST2

This is a re-implementation of Alignment to Median Smoothed Template (AMST), a stack alignment for FIB-SEM datasets, 
originally published here: 

 - [hennies et al. 2020](https://www.nature.com/articles/s41598-020-58736-7)
 - [original github repo](https://github.com/jhennies/amst)

AMST2 uses the same concept compared to the original AMST workflow, while the implementation is now based on Simple 
Elastix (Simple ITK) and includes a pre-alignment workflow.
To reduce the amount of dependencies as well as to increase robustness of the workflow, the use of SIFT was replaced
by Simple Elastix functionality.

## Installation

Clone squirrel and AMST2 repositories:
```shell
cd path/to/src
git clone https://github.com/jhennies/squirrel
git clone https://github.com/jhennies/AMST2
```

Install the packages:
```shell
cd path/to/src
mamba create -n amst2-env -c bioconda -c conda-forge python=3.11 snakemake=8.27
conda activate amst2-env
pip install -e squirrel
pip install -e AMST2
pip install SimpleITK-SimpleElastix
pip install transforms3d
mamba install -c conda-forge vigra
mamba install -c conda-forge opencv
mamba install -c conda-forge zarr=2
```

To use a slurm cluster:
```shell
pip install snakemake-executor-plugin-slurm
```

## Simplified installation

```shell
mamba create -n amst2-env -c bioconda -c conda-forge --override-channels python=3.11 nibabel napari pyqt opencv zarr=2 vigra pandas snakemake=8 snakemake-executor-plugin-slurm
conda activate amst2-env
pip install SimpleITK-SimpleElastix
pip install transforms3d
pip install -e squirrel
pip install -e AMST2
```

## Usage

To test the installation, we recommend using this dataset: https://www.ebi.ac.uk/empiar/EMPIAR-10311/

You can select only the first 32 tif slices (slice_0000.tif to slice_0031.tif) for download, in order to generate a suitable test dataset. 
The examples were tested specifically with this fraction of the dataset, thus this should work properly.
 
Check out the most recent example for [pre alignment and AMST](examples/simple_pre_alignment_and_amst.md)


