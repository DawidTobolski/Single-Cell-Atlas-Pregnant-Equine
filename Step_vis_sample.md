
# Visualization Sample Code: Seurat Object Analysis

This document describes some example visualizations generated during the analysis of single-cell RNA-seq data using a Seurat object. The visualizations are part of a larger analysis process and have been generated for specific markers and clusters.

## Marker Feature Plots

The code iterates over a list of markers (genes) and creates visual plots for each using the `FeaturePlot` function in Seurat. These plots are saved as SVG files. The visualizations include:
- Customized color scaling.
- Removal of axis labels, ticks, and grid lines for a clean presentation.
- Titles with marker names displayed in custom colors.

## Clustering Visualization

The UMAP reduction method is used to visualize clusters of cells based on their transcriptional profiles. The plots show:
- **Cluster IDs**: Different clusters of cells with distinct labels.
- **Group IDs**: Samples split into different conditions (e.g., p33 vs. p42).
The results are saved as SVG files.

## Dot Plot for Gene Expression

A customized dot plot was generated to show the expression of selected genes across clusters. The genes included in the plot are related to immune functions. The resulting dot plot illustrates how gene expression differs across cell types.

## Customized UMAP with Labels

This step generates a UMAP plot with cluster labels and colors assigned to specific cell types. The labels and colors are defined manually to ensure the correct interpretation of the cell types:
- **Epithelial cells** are colored green.
- **Stromal fibroblasts** are colored turquoise.
- **Lymphocytes** are colored lightsalmon.

The customized UMAP plot is saved in SVG format for further analysis.

## Note

These are sample visualizations and part of a larger set of analyses. Only a portion of the visualizations and modifications were transferred to the repository for illustration purposes.
