
# Load necessary libraries for single-cell RNA sequencing data processing and analysis.
library(Seurat)
library(ggplot2)
library(tidyverse)
library(DoubletFinder)
library(harmony)
library(openxlsx)
library(tibble)
library(DESeq2)

# Load Seurat object
p33vsp42.harmony <- readRDS("backup_folder/p33vsp42_harmony.rds")

# Generate UMAP plots for clusters and conditions, saving as SVG files
clusters <- DimPlot(p33vsp42.harmony, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE, raster = FALSE)
condition <- DimPlot(p33vsp42.harmony, reduction = 'umap', group.by = 'group', raster = FALSE)
condition | clusters

# Save the UMAP plots to SVG files
grDevices::svg("clusters.svg")
clusters
dev.off()

grDevices::svg("condition.svg")
condition
dev.off()

# Step 1: Subsetting non-immune clusters
clusters_to_keep <- c(0, 2, 3, 4, 6, 8, 9, 10, 12, 13, 14, 16, 17, 20, 22, 26)

# Subset the non-immune clusters
subset_obj <- subset(p33vsp42.harmony, subset = seurat_clusters %in% clusters_to_keep)

# Process the subset: Find variable features, PCA, and UMAP
subset_obj <- FindVariableFeatures(subset_obj)
subset_obj <- RunPCA(subset_obj)
subset_obj <- RunUMAP(subset_obj, dims = 1:15)

# Perform batch effect correction using Harmony
subset_obj.harmony <- subset_obj %>% RunHarmony(group.by.vars = 'group', plot_convergence = FALSE)

# Check Harmony reductions
subset_obj.harmony@reductions

# Extract Harmony embeddings
subset_obj.harmony.embed <- Embeddings(subset_obj.harmony, "harmony")
print(subset_obj.harmony.embed[1:10, 1:10])

# UMAP, FindNeighbors, and clustering based on Harmony reduction
subset_obj.harmony <- subset_obj.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5)

# Save the non-immune cluster object to an RDS file
saveRDS(subset_obj.harmony, "subset_obj.harmony.non_imune.rds")

# Step 2: Subsetting immune clusters
clusters_to_keep <- c(1, 5, 7, 11, 15, 18, 19, 21, 23, 24, 25, 27)

# Subset the immune clusters
subset_obj <- subset(p33vsp42.harmony, subset = seurat_clusters %in% clusters_to_keep)

# Process the subset: Find variable features, PCA, and UMAP
subset_obj <- FindVariableFeatures(subset_obj)
subset_obj <- RunPCA(subset_obj)
subset_obj <- RunUMAP(subset_obj, dims = 1:15)

# Perform batch effect correction using Harmony
subset_obj.harmony <- subset_obj %>% RunHarmony(group.by.vars = 'group', plot_convergence = FALSE)

# Check Harmony reductions
subset_obj.harmony@reductions

# Extract Harmony embeddings
subset_obj.harmony.embed <- Embeddings(subset_obj.harmony, "harmony")
print(subset_obj.harmony.embed[1:10, 1:10])

# UMAP, FindNeighbors, and clustering based on Harmony reduction
subset_obj.harmony <- subset_obj.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.7)

# Save the immune cluster object to an RDS file
saveRDS(subset_obj.harmony, "subset_obj.harmony.all_imune.rds")
