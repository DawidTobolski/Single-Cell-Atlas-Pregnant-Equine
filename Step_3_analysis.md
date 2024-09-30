
# Step 3: Analysis of Single-Cell RNA-seq Data

## Overview
This script performs single-cell RNA-seq analysis on equine endometrium datasets using the Seurat and Harmony packages. The analysis involves subsetting specific clusters for immune and non-immune cells, followed by batch effect correction, dimensionality reduction, clustering, and visualization.

## Steps:
1. **Load Required Libraries**:
   - Various R libraries are loaded, including Seurat for single-cell analysis, ggplot2 for plotting, and Harmony for batch correction.

2. **Load Seurat Object**:
   - The main Seurat object (`p33vsp42.harmony`) containing UMAP and clustering data is loaded from an `.rds` file.

3. **Generate UMAP Plots**:
   - UMAP plots are generated to visualize the clusters and conditions, which are saved as SVG files.

4. **Subsetting Non-Immune Clusters**:
   - Specific clusters considered non-immune are selected based on predefined cluster indices.
   - The subset undergoes normalization, PCA, and UMAP visualization.
   - Harmony is used for batch effect correction, and the final clusters are saved as an `.rds` file.

5. **Subsetting Immune Clusters**:
   - A second subset of immune clusters is processed similarly with Harmony correction and UMAP visualization.
   - The processed immune cluster data is saved as an `.rds` file.

## Outputs:
- UMAP plots for both clusters and conditions (`clusters.svg`, `condition.svg`).
- Two `.rds` files:
   - `subset_obj.harmony.non_imune.rds`: Non-immune cluster object.
   - `subset_obj.harmony.all_imune.rds`: Immune cluster object.

## Usage
Run the script in an R environment with the necessary dependencies installed. Ensure the Seurat object `p33vsp42.harmony.rds` is located in the `backup_folder` directory.

```
Rscript Step_3_analysis.R
```

This will process the data, generate plots, and save the results into `.rds` files for further analysis.
