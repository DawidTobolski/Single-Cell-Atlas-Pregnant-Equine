
# Differential Gene Expression (DEG) Analysis Between Two Groups

## Overview

This analysis compares gene expression differences between two experimental conditions: **Group A (p33)** and **Group B (p42)**. The data comes from single-cell RNA-seq and has been processed with Seurat and analyzed with DESeq2.

The steps include:
1. **Loading the Processed Data:** Seurat objects containing pre-processed single-cell RNA-seq data were loaded from `.rds` files.
2. **Generating UMAP Visualization:** UMAP was used to visualize clusters, and cluster plots were saved as SVG files.
3. **Differential Gene Expression Analysis (DEG):** Using DESeq2, gene expression differences were compared between two conditions. Both raw counts and normalized counts were provided in the output.

## DEG Results

For each identified cluster, DEG analysis was performed, comparing gene expression levels between **p33** and **p42** conditions. The results include:
- **DE Results:** Table with gene names, log2 fold changes, p-values, and adjusted p-values.
- **Raw Counts:** Unnormalized gene expression counts.
- **Normalized Counts:** Gene expression counts normalized by DESeq2.

Results are saved as `.xlsx` files, with separate sheets for DE results and count data.

## UMAP Clusters Plot

A UMAP plot visualizing cell clusters is saved as `clusters_plot.svg` in the `all_immune` directory.

## DEG Output

All DEG results and count data are saved in the `deg_new_clast` directory. Each cluster has its results saved in a separate Excel file.
