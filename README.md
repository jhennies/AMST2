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

TODO conda recipe

```
conda install -c bioconda snakemake=8
```

## Usage

 1. Conversion to ome.zarr
 2. Template matching alignment
 3. Stack alignment with Simple Elastix
 4. Pre-alignment by combining Template matching (step 3) and Elastix stack alignment (step 4)
 5. Alignment to median smoothed template (AMST)
