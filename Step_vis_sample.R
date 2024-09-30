
# Load necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Step 1: Plot Marker Features
# Create feature plots for a list of markers and save them as SVG files

# Assuming `markers_df` is a data frame containing markers and their associated colors
for (i in 1:nrow(markers_df)) {
  marker <- markers_df$marker[i]
  color <- markers_df$color[i]
  
  # Create FeaturePlot for the marker with customized settings
  p <- FeaturePlot(p33vsp42.harmony, features = marker, cols = c("lightgray", color), 
                   pt.size = 1, min.cutoff = global_min, max.cutoff = global_max) +
       theme_minimal() +
       theme(legend.position = "none",
             axis.title.x = element_blank(),
             axis.title.y = element_blank(),
             axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             plot.background = element_blank(),
             panel.background = element_blank(),
             panel.border = element_blank()) +
       labs(title = marker) +
       theme(plot.title = element_text(color = color, hjust = 0.5, size = 24))
  
  # Save plot as an SVG file
  filename <- paste0("FeaturePlot_", marker, ".svg")
  ggsave(filename, plot = p, width = 6, height = 6, device = "svg")
}

# Step 2: Plot Clusters
# Create UMAP plots for clusters and conditions, and save them as SVG files

clusters <- DimPlot(subset_obj.harmony, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE, raster = FALSE)
condition <- DimPlot(subset_obj.harmony, reduction = 'umap', group.by = 'group', raster = FALSE)

# Save the plots
grDevices::svg("clusters.svg")
print(clusters)
dev.off()

grDevices::svg("condition.svg")
print(condition)
dev.off()

# Step 3: Customized Dot Plot for Specific Genes
# Define a list of genes and generate a dot plot showing their expression

# Create a subset of the Seurat object based on specific clusters
clusters_to_keep <- c("1", "2", "6", "3", "0", "18", "13", "8", "7", "5", "16", "19", "11", "21", "22", "20")
subset_obj.harmony.sub <- subset(subset_obj.harmony, idents = clusters)

# Add information for group and new cluster identifiers
subset_obj.harmony.sub$seurat_clusters <- paste0(subset_obj.harmony.sub$seurat_clusters, "_", subset_obj.harmony.sub$group)
subset_obj.harmony.sub$seurat_clusters <- factor(subset_obj.harmony.sub$seurat_clusters)

# List of specific genes to visualize
new_genes <- c("ENSECAG00000018937", "ENTPD1", "TMIGD2", "GZMA", "CD7", "CD2", "CD3G", "CD3E", "IFNG")

# Create DotPlot and save it as an SVG file
dotplot <- DotPlot(subset_obj.harmony.sub, features = new_genes) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("Immune_DotPlot.svg", plot = dotplot, device = 'svg', width = 20, height = 11)

# Step 4: UMAP with Custom Cluster Labels and Colors
cluster_labels <- c(
  "0" = "0: Epithelial cells", "6" = "6: Epithelial cells", "8" = "8: Epithelial cells",
  "2" = "2: Stromal fibroblasts", "5" = "5: CD3+ lymphocytes"
)

cluster_colors <- c(
  "0" = "mediumseagreen", "6" = "mediumseagreen", "8" = "mediumseagreen",
  "2" = "turquoise", "5" = "lightsalmon"
)

clusters <- DimPlot(p33vsp42.harmony, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE, pt.size = 2, raster = FALSE) +
            scale_color_manual(values = cluster_colors, labels = cluster_labels) +
            guides(color = guide_legend("Cell Types & Clusters", ncol = 1))

# Save the customized UMAP plot
ggsave("clusters_with_custom_legend.svg", plot = clusters, device = "svg", width = 18, height = 14)
