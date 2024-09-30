
# Single-cell RNA-seq Analysis Pipeline for Equine Endometrium Data

## Overview
This repository contains the full workflow and code for processing and analyzing single-cell RNA sequencing (scRNA-seq) data collected from equine endometrium samples. The project involves performing quality control, filtering, normalization, clustering, batch effect correction, and visualization of scRNA-seq data. The results include identifying cell populations and performing differential gene expression analysis.

## Anonymization
To maintain the privacy of research subjects and sensitive data, all sample names and identifiers have been anonymized. Any potential personal or sensitive information has been removed from the data, and the scripts have been adapted accordingly.

## Software Requirements
The following software and R packages are required to run the analysis pipeline:
- **R** (version 4.0 or later)
- `Seurat` for scRNA-seq data processing and analysis
- `DoubletFinder` for detecting doublets in single-cell data
- `harmony` for batch effect correction
- `DESeq2` for differential expression analysis
- `ggplot2`, `tidyverse`, `openxlsx`, `tibble` for visualization and data handling

## Steps in the Analysis

### Step 1: Data Import and Initialization
We begin by importing scRNA-seq data stored in `.h5` format. This data is read and converted into Seurat objects for downstream analysis. The `Read10X_h5` function is used to load the data, and Seurat objects are initialized with minimum filtering parameters.

### Step 2: Quality Control, Filtering, and Normalization
Each dataset undergoes quality control, filtering, and normalization. Cells are filtered based on RNA counts, mitochondrial content, and gene features. Principal Component Analysis (PCA) is performed to reduce dimensionality, followed by clustering and UMAP visualization. DoubletFinder is employed to detect and remove doublets from the data.

### Step 3: Merging and Data Aggregation
Filtered Seurat objects are merged into a combined dataset for further analysis. This allows us to compare and aggregate data across multiple samples, which is critical for downstream analyses such as differential gene expression.

### Step 4: Data Visualization
Various visualization methods are used to evaluate the quality and characteristics of the data:
- Violin plots (`VlnPlot`) for visualizing the distribution of gene counts.
- Feature scatter plots (`FeatureScatter`) for identifying relationships between RNA counts and features.

### Step 6: Saving and Exporting Results
The processed data and metadata are exported to `.xlsx` files using the `writexl` package. Seurat objects are saved as `.rds` files to provide easy access to the processed data for further analyses.

## File Structure
- `analysis.R`: Main script containing the complete workflow.
- `Step_2`: This documentation file.

## Contact Information
For any issues or questions related to this repository, please contact the project team via email at [tobola28@gmail.com].
