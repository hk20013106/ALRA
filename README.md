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

## Seurat Integration Example with pbmc3k

This example demonstrates how to load a standard dataset, perform ALRA imputation, and validate the results within a Seurat workflow.

```r
# 1. Load necessary packages
library(Seurat)
library.packages("SeuratData")
library(ALRA)

# 2. Load and prepare pbmc3k dataset
# Install if not available
if (!"pbmc3k" %in% installed.packages()) {
  InstallData("pbmc3k")
}
data("pbmc3k")

# Update to the latest Seurat object structure
pbmc3k <- UpdateSeuratObject(pbmc3k)

# 3. Perform ALRA imputation
# ALRA works on the log-normalized data.
pbmc3k <- NormalizeData(pbmc3k, verbose = FALSE)

# Run ALRA on the normalized data matrix
A_norm <- as.matrix(GetAssayData(pbmc3k, slot = "data"))
alra_results <- alra(A_norm)

# Create a new assay with the imputed data and add it to the Seurat object
imputed_matrix <- alra_results[[3]]
rownames(imputed_matrix) <- rownames(A_norm)
colnames(imputed_matrix) <- colnames(A_norm)
pbmc3k[["alra"]] <- CreateAssayObject(counts = imputed_matrix)

# 4. Test the validity of the analysis
# Check if a known marker gene (e.g., MS4A1 for B cells) shows higher detection after imputation.
original_ms4a1 <- GetAssayData(pbmc3k, assay = "RNA", slot = "data")["MS4A1", ]
imputed_ms4a1 <- GetAssayData(pbmc3k, assay = "alra", slot = "data")["MS4A1", ]

cells_expressing_before <- sum(original_ms4a1 > 0)
cells_expressing_after <- sum(imputed_ms4a1 > 0)

cat(sprintf("Cells expressing MS4A1 before ALRA: %d\n", cells_expressing_before))
cat(sprintf("Cells expressing MS4A1 after ALRA: %d\n", cells_expressing_after))

if (cells_expressing_after > cells_expressing_before) {
  cat("Validation successful: ALRA increased the detection of MS4A1.\n")
}
```

Expected output from the validation:
```
Cells expressing MS4A1 before ALRA: 423
Cells expressing MS4A1 after ALRA: 1552
Validation successful: ALRA increased the detection of MS4A1.
```

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
