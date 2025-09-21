# Adaptively-thresholded Low Rank Approximation (ALRA)

ALRA provides zero-preserving low-rank completion for single-cell RNA-seq count matrices. The package exposes the core `alra()` pipeline, utilities for rank selection and normalisation, and a Seurat v5-compatible helper (`alraSeurat()`) for seamless integration with modern workflows.

## Key Features
- **Zero-preserving matrix completion** using randomized SVD and adaptive thresholding.
- **Automatic rank heuristics** via `choose_k()` to pick an informative latent dimension.
- **Convenience normalisation** with `normalize_data()` for library-size scaling and log conversion.
- **Seurat integration** through `alraSeurat()` (v5 layers) and the legacy `alraSeurat2()` helper for old slot-based objects.
- **Optional MKL acceleration** by enabling `use.mkl = TRUE` when the [`fastRPCA`](https://github.com/KlugerLab/rpca-mkl) backend is installed.

## Installation

ALRA targets R \>= 4.1 and requires the `Matrix`, `rsvd`, and `SeuratObject` packages. To install the development version from this repository:

```r
# install.packages("remotes")  # if remotes is missing
remotes::install_github("KlugerLab/ALRA", dependencies = TRUE)
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

1. **Prepare a cells-by-genes matrix.** ALRA expects rows = cells, columns = genes.
2. **Normalise the matrix.** Use `normalize_data()` or a custom method.
3. **Select rank `k`.** Call `choose_k()` for an automated heuristic or pass a fixed value.
4. **Run ALRA.** Execute `alra()` and extract the completed matrix (`[[3]]`).

```r
library(ALRA)

# assume counts is a cells x genes matrix
A_norm <- normalize_data(counts)
k_choice <- choose_k(A_norm)
result <- alra(A_norm, k = k_choice$k)
completed_matrix <- result[[3]]
```

`completed_matrix` preserves the original dimensionality while imputing low-abundance transcripts in a zero-aware manner.

## Seurat Integration

ALRA ships with `alraSeurat()` for Seurat v5+ layer objects and `alraSeurat2()` for legacy v2 slot objects.

```r
library(Seurat)
library(ALRA)

pbmc <- NormalizeData(pbmc)
pbmc <- alraSeurat(pbmc, assay = SeuratObject::DefaultAssay(pbmc))
```

Key arguments:
- `assay`: target assay (defaults to `DefaultAssay(object)`).
- `input_layer` / `output_layer`: specify which layer to read/write (`"data"` by default).
- `make_sparse`: keep data sparse (`TRUE` by default).
- Additional parameters are forwarded to `alra()`.

If you load older Seurat v2 objects, call `alraSeurat2()` to update the `@data` slot:

```r
seurat_v2 <- alraSeurat2(seurat_v2)
```

## Worked Example with Included Data

The package bundles a small NK/B-cell dataset (`b_nk_example`) and matching labels (`labels_example`). The following script walks through normalisation, rank selection, completion, and quality checks.

```r
library(ALRA)
library(Matrix)

# 1. Load example data (cells x genes matrix)
data("b_nk_example")
data("labels_example")

# 2. Normalise and pick rank
A_norm <- normalize_data(b_nk_example)
k_choice <- choose_k(A_norm)
cat(sprintf("Chosen rank: %d\n", k_choice$k))

# 3. Run ALRA with the selected rank
alra_out <- alra(A_norm, k = k_choice$k)
A_imputed <- alra_out[[3]]

# 4. Inspect markers recovered by ALRA
marker_summary <- aggregate(
  A_imputed[, c("NCAM1", "CR2")],
  by = list(cell_type = labels_example),
  FUN = function(x) mean(x > 0)
)
print(marker_summary)
```

Expected behaviour:
- `choose_k()` prints the selected rank and intermediate statistics.
- `alra()` reports progress (SVD, thresholding, scaling) and returns a list of three matrices.
- The aggregated marker summary shows improved detection rates for NK marker `NCAM1` and B-cell marker `CR2` relative to the raw data.

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

- **Installation warnings about ORCID**: update the `Authors@R` field in `DESCRIPTION` with a valid ORCID to silence the message.
- **`fastRPCA` missing**: install via GitHub or run with `use.mkl = FALSE` (default) to use the pure R backend.
- **Layer vs. slot errors**: ensure you call `alraSeurat()` for Seurat >= v3/v5 objects. The helper automatically falls back to slot-based APIs when layers are unavailable.

## Citation

If you use ALRA in your research, please cite the accompanying preprint: *Zero-preserving imputation of scRNA-seq data using low-rank approximation* (Linderman et al., bioRxiv 2018).

## License

ALRA is distributed under the MIT license. See `LICENSE.txt` for details.
