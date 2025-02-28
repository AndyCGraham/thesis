# README

This repository contains the code used to produce the thesis of Andrew C Graham. The code is organized into several scripts and notebooks that should be run in a specific order to reproduce the results presented in the thesis.

## Folder Structure

- `scripts/`: Contains scripts for data analysis.
- `data/`: Contains raw and processed data files. (Data will become available upon publication of this work)

## File Run Order

The multiome analysis files were run in the following order:

1. `scripts/multiome/bash/run_process_nf.sh` - Process multiome data from raw reads to RNA and ATAC count matrices for each sample (using cellranger-arc), corrected for ambient RNA contamination with cellbender.
2. `scripts/multiome/R/seuratQC.rmd` - Basic QC and consensus clustering of scMultiome Data
3. `scripts/multiome/python/annotationTRVAE.ipynb` - Annotation of scMultiome Data using Allen Brain atlas reference dataset
4. `scripts/multiome/R/seuratAnalysis.rmd` - In depth analysis of different celltypes
5. `scripts/multiome/R/PRINT_footprinting.rmd` - Footprinting analysis of DG granule cells, required for GRN inference
6. `scripts/multiome/R/CA1_PRINT_footprinting.rmd` - Footprinting analysis of CA1 pyramidal neurons, required for GRN inference
7. `scripts/multiome/R/CA1_GRN_inference.Rmd` - GRN inference of CA1 pyramidal neurons
8. `scripts/multiome/R/DG_GRN_inference.Rmd` - GRN inference of DG granule cells
9. `scripts/multiome/python/SCENIC+.ipynb` - GRN inference with SCENIC+
10. `scripts/multiome/bash/preprocess_signac_for_dictys.sh` - Prepare data for GRN inference with dictys
11. `scripts/multiome/nf/run_dictys_nf.sh` - GRN inference with dictys

```mermaid
    graph TD;
        A[run_process_nf.sh] --> B[seuratQC.rmd];
        B --> C[annotationTRVAE.ipynb];
        C --> D[seuratAnalysis.rmd];
        D --> E[PRINT_footprinting.rmd];
        D --> F[CA1_PRINT_footprinting.rmd];
        D --> G[SCENIC+.ipynb];
        D --> H[preprocess_signac_for_dictys.sh];
        H --> I[run_dictys_nf.sh];
        E --> J[CA1_GRN_inference.Rmd];
        F --> K[DG_GRN_inference.Rmd];
```

Statistical analysis of imaging results is in file `scripts/imaging/imaging_analysis.rmd`.

Statistical analysis of behavioural test results is in file `scripts/behaviour/behaviour_analysis.rmd`.

## Data Availability

The data used in this analysis will become available upon the publication of this work. 

For any questions or issues, please contact Andrew C Graham.
