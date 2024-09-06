# SCCAF-D: Single-Cell Analysis Framework for Deconvolution
**SCCAF-D** is a novel method for **cell type deconvolution** in single-cell RNA-seq data. By integrating multiple single-cell datasets and employing advanced machine learning techniques, SCCAF-D produces optimized reference data to ensure reliable and accurate deconvolution results. This approach provides valuable insights into the cellular composition and heterogeneity of biological samples.

## Metadata Requirements

To perform deconvolution, the metadata of the single-cell datasets should contain the following essential columns:

- **cellID**: Unique identifiers for individual cells in the dataset.
- **cellType**: Annotations or labels indicating the type of each cell.
- **sampleID**: Identifiers for biological samples to facilitate comparative analysis across different conditions.

These three columns are critical for accurate cell type deconvolution using SCCAF-D.

## Parameters

SCCAF-D requires the following parameters to perform deconvolution:
   - **Bulk**: The bulk RNA-seq data used to be deconvolved.
   - **Reference**: The reference data used for deconvolution.
   - **Transformation** (string): The method for transforming both the bulk and reference data. Options include `none` (default), `log`, `sqrt`, and `vst`.
   - **Deconv_type** (string): Deconvolution methods are categorized into bulk and single-cell (`sc`) approaches based on the reference source. Bulk methods, such as `CIBERSORT` and `FARDEEP`, use a reference signature from sorted cell types or a marker gene list. In contrast, single-cell (`sc`) methods, such as `MuSiC` and `DWLS`, use single-cell datasets as reference data.
   - **Normalization_C** (string): Normalization for reference data. Eighteen normalization methods are supported, including `column`, `row`, `mean`, `column z-score`, `global z-score`, `column min-max`, `global min-max`, `LogNormalize`, `none`, `QN`, `TMM`, `UQ`, `median ratios`, `TPM`, `SCTransform`, `scran`, `scater`, and `Linnorm`.
   - **Normalization_T** (string): Normalization for bulk data. The same normalization methods as listed for reference data. **Suggestion**: If using the single-cell (`sc`) method for deconvolution, the choice for this parameter should match the reference data normalization.
   - **Method** (string): Twenty-five deconvolution algorithms are available, including `DWLS`, `FARDEEP`, `MuSiC`, `nnls`, `RLR`, `EpiDISH`, `OLS`, `EPIC`, `elasticNet`, `lasso`, `proportionsInAdmixture`, `ridge`, `CIBERSORT`, `SCDC`, `BisqueRNA`, `CDSeq`, `CPM`, `DCQ`, `DSA`, `DeconRNASeq`, `TIMER`, `deconf`, `dtangle`, `ssFrobenius`, and `ssKL`.
   - **Number_cells** (integer): The number of cells to select when preparing simulated 'pseudobulk' from single-cell data.
   - **To_remove** (string): Specify any cell types to exclude from the reference data.
   - **Num_cores** (integer): The number of cores to use for parallelization during deconvolution.
   - **NormTrans** (logical): Whether to perform data normalization or transformation first. `T` indicates normalization first, while `F` indicates transformation first.
   - **Return_expr** (logical): Whether to return the estimated expression matrix of the bulk data. Default is `FALSE`.
   - **Batch_key** (string): The parameter used to calculate highly variable genes in `SCANPY`.
   - **Span** (numeric): The fraction of cells used when estimating variance in the loess model fit in `SCANPY` (when `flavor = 'seurat_v3'`).
   - **Python_home** (string): The path to the Python executable.

----

## Installation

To install the required SCCAF package, use one of the following commands:

```shell
conda install sccaf
```

or

```shell
pip install sccaf
```

Additionally, several R and Bioconductor packages are needed:

```R
# CRAN packages
packages <- c("devtools", "BiocManager", "data.table", "ggplot2", "tidyverse", 
              "reticulate", "pheatmap", "Matrix", "matrixStats", "gtools",
              "foreach", "doMC", "doSNOW", "Seurat", "sctransform", "nnls", 
              "MASS", "glmnet")
install.packages(packages)

# Bioconductor packages
bioc_packages <- c('limma', 'edgeR', 'DESeq2', 'pcaMethods', 'BiocParallel', 
                   'preprocessCore', 'scater', 'SingleCellExperiment', 'Linnorm',
                   'DeconRNASeq', 'multtest', 'GSEABase', 'annotate', 'genefilter', 
                   'graph', 'MAST', 'Biobase', 'sparseMatrixStats')
BiocManager::install(bioc_packages)
```

## Example Usage in R

```R
# Set working directory to the location of your data
setwd('***')

# Load the SCCAF-D R function
source('./SCCAF_D.R')

# Define the parameters for deconvolution
param <- c("bulk.rds", "single-reference.rds", "none", "sc", "TMM", "TMM", "DWLS", 
           10000, "none", 1, 'T')

# Specify the path to the Python environment
python_home <- '/home/feng_shuo/miniconda3/envs/sccaf/bin/python'

# Run SCCAF-D deconvolution
results <- SCCAF_D(param, python_home = python_home)
```

## Citation
To cite SCCAF-D, please refer to the following:

> Feng *et al.* "Alleviating batch effects in cell type deconvolution." (Preprint available at Research Square)