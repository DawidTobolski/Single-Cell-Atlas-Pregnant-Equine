
# Load necessary libraries for single-cell RNA sequencing data processing and analysis.
library(Seurat)          # For single-cell RNA-seq data analysis.
library(ggplot2)         # For data visualization.
library(tidyverse)       # A collection of R packages for data manipulation and visualization.
library(DoubletFinder)   # For detecting doublets in single-cell RNA sequencing data.
library(harmony)         # For batch effect correction in single-cell RNA sequencing data.
library(openxlsx)        # For exporting results into Excel format.
library(tibble)          # For better data handling in tibbles (modern data frames).
library(DESeq2)          # For differential gene expression analysis.

# Set a random seed for reproducibility in data analysis.
set.seed(1234)

# Step 1: Data Import (Anonymized File Names)
# Anonymized list of HDF5 files containing the single-cell RNA sequencing data.
# The files correspond to different experimental samples. Each file name has been anonymized to 
# remove identifying information.
file_list <- c("s01_13P42.h5", "s02_2P33.h5", "s02_14P42.h5", "s4_4P33.h5", "s06_6P33.h5", "s12_12P42.h5")

# Loop through the list of anonymized HDF5 files and perform the following:
# - Read the data from the HDF5 files using Read10X_h5.
# - Extract the project name from each file by removing the '.h5' suffix.
# - Create Seurat objects for each dataset, which store the single-cell data for analysis.
for (file_name in file_list) {
  data <- Read10X_h5(filename = file_name)   # Read the HDF5 file.
  project_name <- gsub(".h5", "", file_name) # Remove '.h5' extension to get project name.
  
  # Store the data for each project in variables, appending ".data" to the project name.
  assign(paste0(project_name, ".data"), data)
  
  # Create a Seurat object for each dataset.
  obj_name <- paste0(project_name, ".obj")
  counts <- get(paste0(project_name, ".data"))
  obj <- CreateSeuratObject(counts = counts, project = project_name, min.cells = 3, min.features = 200)
  
  # Store the Seurat object for further analysis.
  assign(obj_name, obj)
}

# Step 2: Quality Control, Filtering, and Normalization
# This loop performs quality control, filtering, normalization, and dimensionality reduction 
# for each Seurat object created in Step 1.

for (file_name in file_list) {
  project_name <- gsub(".h5", "", file_name)  # Extract the project name by removing the '.h5' extension.
  obj_name <- paste0(project_name, ".obj")    # Name for the original Seurat object.
  obj_name.filtered <- paste0(project_name, ".filtered")  # Name for the filtered Seurat object.
  
  filtered <- get(obj_name)  # Retrieve the Seurat object.
  
  # Add percentage of mitochondrial and ribosomal RNA reads to metadata.
  # Mitochondrial RNA is often a marker of low-quality cells.
  filtered$mitoPercent <- PercentageFeatureSet(filtered, pattern = '^MT-')  # Mitochondrial genes.
  filtered$riboPercent <- PercentageFeatureSet(filtered, pattern = '^Rn')   # Ribosomal genes.
  
  # Filter cells based on several quality metrics:
  # - RNA counts (nCount_RNA) between 750 and 15000.
  # - Number of detected genes (nFeature_RNA) between 250 and 6500.
  # - Mitochondrial content (mitoPercent) less than 5%.
  filtered <- subset(filtered, subset = nCount_RNA > 750 & nCount_RNA < 15000 & nFeature_RNA > 250 & nFeature_RNA < 6500 & mitoPercent < 5)
  
  # Normalize the data to account for differences in sequencing depth across cells.
  filtered <- NormalizeData(object = filtered)
  
  # Identify the most variable genes in the dataset for downstream analysis.
  filtered <- FindVariableFeatures(object = filtered)
  
  # Scale the data to ensure that features are comparable.
  filtered <- ScaleData(object = filtered)
  
  # Perform Principal Component Analysis (PCA) to reduce dimensionality.
  filtered <- RunPCA(object = filtered)
  
  # Plot an elbow plot to help identify the optimal number of principal components to use.
  ElbowPlot(filtered)
  
  # Find cell neighbors based on the PCA results, which are used for clustering.
  filtered <- FindNeighbors(object = filtered, dims = 1:20)
  
  # Perform clustering using the Louvain algorithm.
  filtered <- FindClusters(object = filtered)
  
  # Generate UMAP projections for visualization of clusters.
  filtered <- RunUMAP(object = filtered, dims = 1:20)
  
  # Doublet detection using DoubletFinder to identify potential doublet cells.
  sweep.res.list <- paramSweep_v3(filtered, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  
  # Find the optimal value of the pK parameter for DoubletFinder.
  bcmvn <- find.pK(sweep.stats)
  ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
    geom_point() +
    geom_line()
  
  pK <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% select(pK)  # Select optimal pK value.
  pK <- as.numeric(as.character(pK[[1]]))  # Convert pK to a numeric value.
  
  # Homotypic proportion estimation (helps refine doublet prediction).
  annotations <- filtered@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  
  # Estimate the expected number of doublets.
  nExp_poi <- round(0.046 * nrow(filtered@meta.data))  # Based on 4.6% doublet rate.
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))  # Adjust for homotypic proportion.
  
  # Perform doublet detection using DoubletFinder.
  filtered <- doubletFinder_v3(filtered, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  
  # Save the filtered object for further analysis.
  assign(obj_name.filtered, filtered)
}

# Step 3: Seurat Objects and Merging Data
seurat_objects <- list(
  s01_13P42 = s01_13P42.filtered, 
  s02_2P33 = s02_2P33.filtered, 
  s02_14P42 = s02_14P42.filtered, 
  s4_4P33 = s4_4P33.filtered, 
  s06_6P33 = s06_6P33.filtered, 
  s12_12P42 = s12_12P42.filtered
)

# Summary table of features and samples per object
seurat_table.f <- data.frame(Name = character(), Features = numeric(), Samples = numeric(), stringsAsFactors = FALSE)
for (name in names(seurat_objects)) {
  obj <- seurat_objects[[name]]
  features <- obj@assays$RNA@counts %>% ncol()
  samples <- obj@assays$RNA@counts %>% nrow()
  seurat_table.f <- rbind(seurat_table.f, data.frame(Name = name, Features = features, Samples = samples))
}
print(seurat_table.f)

# Step 4: Merging and Preparing Combined Data
p33vsp42 <- merge(
  s01_13P42.single, 
  y = c(s02_2P33.single, s02_14P42.single, s4_4P33.single, s06_6P33.single, s12_12P42.single), 
  add.cell.ids = c("s01_13P42", "s02_2P33", "s02_14P42", "s4_4P33", "s06_6P33", "s12_12P42"), 
  project = "p33vsp42"
)

p33vsp42$group <- plyr::mapvalues(
  x = p33vsp42$orig.ident, 
  from = c("s01_13P42", "s02_2P33", "s02_14P42", "s4_4P33", "s06_6P33", "s12_12P42"), 
  to = c("p42", "p33", "p42", "p33", "p33", "p42")
)

# Step 5: Data Visualization and Summary Plots
VlnPlot(p33vsp42, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(p33vsp42, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')

# Save Metadata to Excel
metadata_df <- p33vsp42@meta.data %>% select(orig.ident, nFeature_RNA)
writexl::write_xlsx(metadata_df, "p33vsp42_metadata.xlsx")

perform_DE_analysis()

# Step 7: Save RDS Object for Backup
saveRDS(p33vsp42.harmony, "backup_folder/p33vsp42_harmony.rds")
