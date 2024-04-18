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
```
cd path/to/src
git clone https://github.com/jhennies/squirrel
git clone https://github.com/jhennies/AMST2
```

Install the packages:
```
cd path/to/src
mamba create -n amst2-env -c conda-forge -c bioconda python=3.11 snakemake=8
conda activate amst2-env
pip install -e squirrel
pip install -e AMST2
pip install SimpleITK-SimpleElastix
pip install transforms3d
mamba install -c conda-forge opencv
mamba install -c conda-forge zarr
```

To use a slurm cluster:
```
pip install snakemake-executor-plugin-slurm
```

## Usage

 1. Conversion to ome.zarr
 2. Template matching alignment
 3. Stack alignment with Simple Elastix
 4. Pre-alignment by combining Template matching (step 3) and Elastix stack alignment (step 4)
 5. Alignment to median smoothed template (AMST)
