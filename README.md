
# Single-Cell Atlas of the Pregnant Equine Endometrium Before and After Implantation

This repository contains the data and analysis code for the study: **"Single-cell atlas of the pregnant equine endometrium before and after implantation"** by Jaworska et al. The study explores the transcriptome and cellular composition of the equine endometrium during early pregnancy, specifically before and after embryo implantation.

## Key Findings
- Endometrial NK cells are the most abundant leukocyte population during implantation in horses.
- Antigen-presenting cells (APCs) are the most responsive immune cell population across implantation.
- Key contributors to the endometrial immune landscape include innate lymphoid cells (ILCs) and mucosal-associated invariant T cells (MAIT cells), identified for the first time in the horse.
- Despite differences in placentation, horses share similarities with humans in the mechanisms of embryo implantation, such as the increase in NK cells and upregulation of CXCR4 expression.

## Analysis Pipeline
The analysis pipeline consists of several steps, each implemented in R scripts:

### Step 1: Data Preparation with Cell Ranger
- Creation of a custom reference genome.
- Processing of raw sequencing data (FASTQ files) to generate count matrices.
- Details and instructions are provided in Step_1_Cell_Ranger.md.

### Step 2: Quality Control and Initial Analysis
- Data import, quality control, filtering, normalization, and dimensionality reduction using Seurat.
- Doublet detection and removal using DoubletFinder.
- Batch effect correction using Harmony.
- Initial clustering and visualization.
- The code is available in Step_2_analysis.R, and a description is in Step_2_analysis.md.

### Step 3: Subsetting and Further Analysis
- Subsetting of immune and non-immune cell clusters.
- Further batch correction, dimensionality reduction, clustering, and visualization of the subsets.
- The code is available in Step_3_analysis.R, and a description is in Step_3_analysis.md.

### Differential Gene Expression Analysis
- Identification of differentially expressed genes (DEGs) between day 33 and day 42 of pregnancy using DESeq2.
- The code is available in Step_deg.R, and a description is in Step_deg.md.

### Visualization
- Generation of various plots, including feature plots, UMAP plots, and dot plots, to visualize gene expression and cell clusters.
- The code is available in Step_vis_sample.R, and a description is in Step_vis_sample.md.

## Data Availability
Sequencing data were uploaded to... [Please fill in the appropriate repository or accession number].

## Software Requirements
- R (version 4.0 or later)
- R packages: Seurat, DoubletFinder, harmony, DESeq2, ggplot2, tidyverse, openxlsx, tibble

## Contact
For any questions or issues, please contact Dawid Tobolski at [tobola28@gmail.com].

## Note
The code provided here has undergone multiple iterations and adjustments to optimize visualizations. Thus, the version in this repository may not represent the exact final version used in the study but is reflective of the core analysis steps.
