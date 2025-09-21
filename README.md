# Adaptively-thresholded Low Rank Approximation (ALRA)

ALRA provides zero-preserving low-rank completion for single-cell RNA-seq count matrices. The package exposes the core `alra()` pipeline, utilities for rank selection and normalisation, and seamless integration with modern Seurat workflows.

## Key Features
- **Zero-preserving matrix completion** using randomized SVD and adaptive thresholding.
- **Automatic rank heuristics** via `choose_k()` to pick an informative latent dimension.
- **Convenience normalisation** with `normalize_data()` for library-size scaling and log conversion.
- **Optional MKL acceleration** by enabling `use.mkl = TRUE` when the [`fastRPCA`](https://github.com/KlugerLab/rpca-mkl) backend is installed.

## Installation

ALRA targets R \>= 4.1 and requires the `Matrix`, `rsvd`, and `SeuratObject` packages. To install the development version from this repository:

```r
# install.packages("remotes")  # if remotes is missing
remotes::install_github("hk20013106/ALRA", dependencies = TRUE)
```

For Windows users leveraging MKL acceleration, install `fastRPCA` before calling `alra()` with `use.mkl = TRUE`:

```r
remotes::install_github("KlugerLab/rpca-mkl/fastRPCA")
```

After installation, load the package as usual:

```r
library(ALRA)
```

## Quick Start: Core Workflow

1. **Prepare a genes-by-cells matrix.** ALRA expects rows = genes, columns = cells.
2. **Normalise the matrix.** Use `Seurat::NormalizeData` or a custom method.
3. **Select rank `k`.** Call `choose_k()` for an automated heuristic or pass a fixed value.
4. **Run ALRA.** Execute `alra()` and extract the completed matrix (`[[3]]`).

```r
library(ALRA)
library(Seurat)

# Assume pbmc is a Seurat object with raw counts
pbmc <- NormalizeData(pbmc, verbose = FALSE)
A_norm <- as.matrix(GetAssayData(pbmc, slot = "data"))

k_choice <- choose_k(A_norm)
result <- alra(A_norm, k = k_choice$k)
completed_matrix <- result[[3]]
```

`completed_matrix` preserves the original dimensionality while imputing low-abundance transcripts in a zero-aware manner.

## Examples

### Example 1: Basic Usage with B and NK cell data

This example demonstrates the core functionality of ALRA on a small, built-in dataset of B and NK cells. We load the data, normalize it, choose a rank `k`, and then run ALRA. Finally, we compare the percentage of non-zero cells for marker genes `NCAM1` (for NK cells) and `CR2` (for B cells) before and after imputation.

The complete, runnable code for this analysis is available in the `b_nk_example.R` script.

**Code:**
```r
# b_nk_example.R

# Load the ALRA package
# Assuming it's installed, otherwise: remotes::install_github("hk20013106/ALRA")
library(ALRA)

# Load the example data included with the package
data("b_nk_example")
data("labels_example")

# Normalize the data
A_norm <- normalize_data(b_nk_example)

# Choose the rank k
k_choice <- choose_k(A_norm)

# Run ALRA
A_norm_completed <- alra(A_norm, k=k_choice$k)[[3]]

# Print statistics before ALRA
print("Before ALRA:")
print(aggregate(A_norm[,c('NCAM1','CR2')], by=list(" "=labels_example),FUN=function(x) round(c(percent=100*sum(x>0)/length(x)),1)))

# Print statistics after ALRA
print("After ALRA:")
print(aggregate(A_norm_completed[,c('NCAM1','CR2')], by=list(" "=labels_example),FUN=function(x) round(c(percent=100*sum(x>0)/length(x)),1)))
```

**Output:**
```
[1] "Before ALRA:"
          NCAM1 CR2
1 b_cells   0.0 0.9
2 cd56_nk   3.9 0.0
[1] "After ALRA:"
          NCAM1  CR2
1 b_cells   0.1 70.4
2 cd56_nk  98.8  0.5
```
As shown in the output, ALRA significantly increases the expression of the correct marker in each cell type, clarifying their identity.

### Example 2: Seurat Integration with pbmc3k

This example demonstrates how to integrate ALRA into a standard Seurat workflow. We use the `pbmc3k` dataset, perform a standard analysis (normalization, scaling, PCA, UMAP) on the original data, and then use ALRA to impute the expression. The key point is that we use the **same UMAP** generated from the original data to visualize both the original and imputed expression, allowing for a direct comparison.

The complete, runnable code for this analysis is available in the `full_analysis_example.R` script.

**Code:**
```r
# full_analysis_example.R (Simplified)

# 1. Load necessary packages
# ... (installation checks omitted for brevity)
library(Seurat)
library(SeuratData)
library(ALRA)
library(patchwork)
library(ggplot2)

# 2. Load and prepare pbmc3k dataset
data("pbmc3k")
pbmc3k <- UpdateSeuratObject(pbmc3k)

# 3. Perform ALRA imputation
pbmc3k <- NormalizeData(pbmc3k, verbose = FALSE)
A_norm <- as.matrix(GetAssayData(pbmc3k, layer = "data"))
alra_results <- alra(A_norm)
imputed_matrix <- alra_results[[3]]
rownames(imputed_matrix) <- rownames(A_norm)
colnames(imputed_matrix) <- colnames(A_norm)
pbmc3k[["alra"]] <- CreateAssayObject(counts = imputed_matrix)

# 4. Standard analysis on original data (to generate one UMAP)
DefaultAssay(pbmc3k) <- "RNA"
pbmc3k <- FindVariableFeatures(pbmc3k, verbose = FALSE)
pbmc3k <- ScaleData(pbmc3k, verbose = FALSE)
pbmc3k <- RunPCA(pbmc3k, verbose = FALSE)
pbmc3k <- FindNeighbors(pbmc3k, dims = 1:10, verbose = FALSE)
pbmc3k <- FindClusters(pbmc3k, resolution = 0.5, verbose = FALSE)
pbmc3k <- RunUMAP(pbmc3k, dims = 1:10, verbose = FALSE, reduction.name = "umap.rna")

# 5. Generate and save side-by-side comparison plots
markers_to_plot <- c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")
plot_pairs <- list()

for (gene in markers_to_plot) {
    p_rna <- FeaturePlot(pbmc3k, features = gene, reduction = "umap.rna", order = TRUE) + NoLegend() + NoAxes() + ggtitle(paste(gene, "(Original)"))
    
    # Create a temporary Seurat object for ALRA visualization to avoid altering the main object's state
    temp_seurat <- pbmc3k
    DefaultAssay(temp_seurat) <- "alra"
    
    p_alra <- FeaturePlot(temp_seurat, features = gene, reduction = "umap.rna", order = TRUE) + NoAxes() + ggtitle(paste(gene, "(ALRA)"))
    
    plot_pairs[[gene]] <- p_rna | p_alra
}

combined_figure <- wrap_plots(plot_pairs, ncol = 2)
ggsave("alra_marker_comparison.png", plot = combined_figure, width = 10, height = 20, dpi = 300, limitsize = FALSE)
```

**Visual Comparison:**

The figure below shows `FeaturePlots` for several marker genes. Both original and ALRA-imputed expression values are visualized on the *same* UMAP projection, which was calculated from the original data. This clearly demonstrates how ALRA enhances the gene expression signals without altering the underlying cellular structure.

![Comparison of Original and ALRA Imputed Data on the same UMAP](alra_marker_comparison.png)

*Figure: Side-by-side FeaturePlots for canonical marker genes. Both plots in each pair use the same UMAP coordinates derived from the original RNA data. The left plot shows the original expression, and the right plot shows the ALRA-imputed expression. ALRA clarifies expression patterns, making cell populations more distinct.*


## Testing and Reproducibility

Developers can run the included unit tests (requires `pkgload`, `testthat`, `Seurat`, and `SeuratObject`):

```r
pkgload::load_all(export_all = FALSE, helpers = FALSE, quiet = TRUE)
testthat::test_dir("tests/testthat", reporter = "summary")
```

Set a random seed before calling `alra()` to mirror published results, as the randomized SVD initialisation introduces stochasticity:

```r
set.seed(42)
```

## Troubleshooting

- **`fastRPCA` missing**: install via GitHub or run with `use.mkl = FALSE` (default) to use the pure R backend.
- **Layer vs. slot errors**: Ensure your Seurat objects are updated to the v5 structure using `UpdateSeuratObject()`.

## Citation

If you use ALRA in your research, please cite the accompanying preprint: *Zero-preserving imputation of scRNA-seq data using low-rank approximation* (Linderman et al., bioRxiv 2018).

## License

ALRA is distributed under the MIT license. See `LICENSE.txt` for details.
