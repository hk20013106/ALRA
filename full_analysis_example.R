#
# This script demonstrates a full analysis pipeline comparing original data
# with ALRA-imputed data using the pbmc3k dataset from SeuratData.
#

# 1. Load necessary packages
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("SeuratData", quietly = TRUE)) install.packages("SeuratData")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("ALRA", quietly = TRUE)) remotes::install_github("hk20013106/ALRA")

library(Seurat)
library(SeuratData)
library(ALRA)
library(patchwork)
library(ggplot2)

# 2. Load and prepare pbmc3k dataset
if (!"pbmc3k" %in% installed.packages()) {
  InstallData("pbmc3k")
}
data("pbmc3k")
pbmc3k <- UpdateSeuratObject(pbmc3k)

# 3. Perform ALRA imputation
pbmc3k <- NormalizeData(pbmc3k, verbose = FALSE)
A_norm <- as.matrix(GetAssayData(pbmc3k, layer = "data"))

# Run ALRA and create a new assay
alra_results <- alra(A_norm)
imputed_matrix <- alra_results[[3]]
rownames(imputed_matrix) <- rownames(A_norm)
colnames(imputed_matrix) <- colnames(A_norm)
pbmc3k[["alra"]] <- CreateAssayObject(counts = imputed_matrix)

cat("ALRA imputation complete.\n")

# 4. Standard analysis on original data
DefaultAssay(pbmc3k) <- "RNA"
pbmc3k <- FindVariableFeatures(pbmc3k, verbose = FALSE)
pbmc3k <- ScaleData(pbmc3k, verbose = FALSE)
pbmc3k <- RunPCA(pbmc3k, verbose = FALSE)
pbmc3k <- FindNeighbors(pbmc3k, dims = 1:10, verbose = FALSE)
pbmc3k <- FindClusters(pbmc3k, resolution = 0.5, verbose = FALSE)
pbmc3k <- RunUMAP(pbmc3k, dims = 1:10, verbose = FALSE, reduction.name = "umap.rna")

cat("Analysis on original data complete.\n")

# 5. Standard analysis on ALRA-imputed data
DefaultAssay(pbmc3k) <- "alra"
# The data is already imputed and normalized. We scale all genes.
all.genes <- rownames(pbmc3k)
pbmc3k <- ScaleData(pbmc3k, features = all.genes, verbose = FALSE)
pbmc3k <- RunPCA(pbmc3k, features = all.genes, verbose = FALSE, reduction.name = "pca.alra")
pbmc3k <- FindNeighbors(pbmc3k, dims = 1:10, reduction = "pca.alra", verbose = FALSE)
pbmc3k <- FindClusters(pbmc3k, resolution = 0.5, verbose = FALSE, cluster.name = "alra_clusters")
pbmc3k <- RunUMAP(pbmc3k, reduction = "pca.alra", dims = 1:10, verbose = FALSE, reduction.name = "umap.alra")

cat("Analysis on ALRA-imputed data complete.\n")

# 6. Generate and save comparison plots
markers_to_plot <- c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")

# Plot RNA
p_rna <- FeaturePlot(
  pbmc3k,
  features = markers_to_plot,
  reduction = "umap.rna",
  order = TRUE,
  combine = FALSE
)
p_rna <- lapply(p_rna, function(p) p + NoLegend() + NoAxes())
plot_rna_wrapped <- wrap_plots(p_rna, ncol = 3) + plot_annotation(title = "Original Data (RNA Assay)")

# Plot ALRA
p_alra <- FeaturePlot(
  pbmc3k,
  features = markers_to_plot,
  reduction = "umap.alra",
  order = TRUE,
  combine = FALSE
)
p_alra <- lapply(p_alra, function(p) p + NoLegend() + NoAxes())
plot_alra_wrapped <- wrap_plots(p_alra, ncol = 3) + plot_annotation(title = "Imputed Data (ALRA Assay)")

# Combine and save
combined_figure <- plot_rna_wrapped / plot_alra_wrapped
ggsave("alra_marker_comparison.png", plot = combined_figure, width = 12, height = 8, dpi = 300)

cat("Comparison plot saved to alra_marker_comparison.png\n")
