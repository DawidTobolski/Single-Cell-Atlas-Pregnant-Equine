
# Load necessary libraries for analysis
library(Seurat)       # For single-cell RNA-seq analysis
library(DESeq2)       # For differential gene expression analysis
library(openxlsx)     # For exporting results to Excel format

# Set working directory (anonymized path)
setwd("D:/anonymized_path/Analiza_sc_21.11.2023")

# Step 1: Load the pre-processed Seurat object
subset_obj.harmony <- readRDS("anonymized_backup/p33vsp42_harmony.rds")

# Step 2: Generate and save UMAP clusters plot
clusters <- DimPlot(subset_obj.harmony, reduction = 'umap', group.by = "ident", label = TRUE, raster = FALSE)

# Create directory if it doesn't exist
if (!dir.exists("all_immune")) {
  dir.create("all_immune")
}

# Save the clusters plot to SVG format
grDevices::svg("all_immune/clusters_plot.svg")
print(clusters)  # Print is necessary for saving the plot
dev.off()

# Step 3: DEG analysis using raw counts
# Set the identity of clusters for DEG analysis
subset_obj.harmony <- SetIdent(subset_obj.harmony, value = "seurat_clusters")

# Update the 'samples' column with new clustering identity
subset_obj.harmony$samples <- paste0(subset_obj.harmony$group, subset_obj.harmony$orig.ident)

# Pseudobulk expression matrix creation based on new clustering
pseudo <- AggregateExpression(subset_obj.harmony, 
                              group.by = c("seurat_clusters", "samples"),
                              assays = 'RNA',
                              slot = "counts",
                              return.seurat = FALSE)

# Transpose and modify the pseudobulk data for analysis
pseudo <- pseudo$RNA
pseudo.t <- t(pseudo)
pseudo.t <- as.data.frame(pseudo.t)

splitRows <- gsub('_.*', '', rownames(pseudo.t))
pseudo.split <- split.data.frame(pseudo.t, f = factor(splitRows))

# Modify row names to prepare data for DE analysis
pseudo.split.modified <- lapply(pseudo.split, function(x){
  rownames(x) <- gsub("^[^_]*_", "", rownames(x))
  t(x)
})

# Function to perform differential expression (DE) analysis
perform_DE_analysis <- function() {
  
  # Create DEG output directory if not already existing
  if (!dir.exists("deg_new_clast")) {
    dir.create("deg_new_clast")
  }
  
  cluster_list <- names(pseudo.split.modified)
  output_dir <- "deg_new_clast"  # Output directory for DEG results
  
  for (cluster in cluster_list) {
    tryCatch({
      # Prepare count matrix for DESeq2
      counts_cluster <- pseudo.split.modified[[cluster]]
      colData <- data.frame(samples = colnames(counts_cluster))
      colData <- colData %>%
        mutate(condition = ifelse(grepl('p33', samples), 'p33', 'p42')) %>%
        column_to_rownames(var = 'samples')
      
      dds <- DESeqDataSetFromMatrix(countData = counts_cluster,
                                    colData = colData,
                                    design = ~ condition)
      
      # Filter genes with low counts
      keep <- rowSums(counts(dds)) >= 100
      dds <- dds[keep,]
      
      # Perform DESeq2 analysis
      dds <- DESeq(dds)
      res <- results(dds, name = "condition_p42_vs_p33", pAdjustMethod = "none")
      res <- cbind(gene = rownames(res), res)
      
      # Add raw and normalized count data for each gene
      raw_counts <- counts(dds, normalized = FALSE)
      raw_counts_df <- as.data.frame(raw_counts)
      raw_counts_df$gene <- rownames(raw_counts_df)
      
      norm_counts <- counts(dds, normalized = TRUE)
      norm_counts_df <- as.data.frame(norm_counts)
      norm_counts_df$gene <- rownames(norm_counts_df)
      
      # Save DE results and raw counts to an Excel file
      write.xlsx(list(DE_Results = res, Raw_Counts = raw_counts_df, Norm_Counts = norm_counts_df),
                 file = paste0(output_dir, "/res_cluster_", cluster, ".xlsx"),
                 rowNames = TRUE)
      
    }, error = function(e) {
      message("An error occurred in cluster ", cluster, ": ", conditionMessage(e))
    })
  }
}

# Run the DEG analysis function
perform_DE_analysis()
