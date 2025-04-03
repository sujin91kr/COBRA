

# COBRA
**A New Batch Effect Correction Pipeline for Single-Cell RNA-Seq Data**

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Introduction
**COBRA** is an R package that provides a novel pipeline for correcting batch effects in single-cell RNA-seq data. It uses an orthogonal projection approach to remove unwanted variation while optionally controlling for known covariates. If cell type labels are unknown, COBRA can estimate pseudo-cells by clustering and iteratively removing batch effects.

### Key Features
- **Batch Effect Removal**: Orthogonal projection-based correction
- **Pseudo-Cell Estimation**: Automatic cell clustering if cell types are not provided
- **Flexibility**: Optionally stratify by cell type or include additional covariates
- **Memory Efficiency**: `reduce` option to handle large datasets

---

## Installation

### 1) Install from GitHub
If you have **devtools** installed:

```r
# install.packages("devtools")
devtools::install_github("sujin91kr/COBRA")
library(COBRA)
```

### 2) Install Locally
1. Clone or download this repository:
   ```bash
   git clone https://github.com/sujin91kr/COBRA.git
   ```
2. In R:
   ```r
   # install.packages("devtools")
   devtools::install("sujin91kr/COBRA")
   library(COBRA)
   ```

---

## Usage

Below is a quick example demonstrating the main function **`BatchCorrect`**.  
We assume you have a count matrix (`data.mat`) with rows = cells and columns = genes, and a metadata data.frame (`meta.mat`) containing a batch label column (`batch`) and optionally a cell type column (`CellType`).

```r
# Example: Basic usage
corrected <- BatchCorrect(
  data.mat    = data.mat,
  meta.mat    = meta.mat,
  batch.var   = "batch",
  cell.var    = "CellType",
  coi.var     = c("Gender", "Age"),  # optional covariates
  stratify    = FALSE,
  reduce      = FALSE,
  sparse.tol  = 1e-6
)
```

**Arguments**:
- `data.mat`: Matrix of expression data (cells x genes)
- `meta.mat`: Data frame with batch labels and optional cell type / covariates
- `batch.var`: Column name in `meta.mat` for batch label
- `cell.var`: Column name for cell type label (if known)
- `coi.var`: A vector of additional covariate names in `meta.mat`
- `stratify`: Logical, whether to stratify batch effect correction by cell type
- `reduce`: Logical, if TRUE uses alternative data orientation for memory efficiency
- `sparse.tol`: Numeric tolerance for near-zero cutoff

**Returned Value**:
- If `cell.var` is given: returns an adjusted expression matrix
- If `cell.var` is `NULL`: returns a list with `(adj, pseudocell)`, where `pseudocell` is the cluster-based pseudo-cell assignment

---

## Additional Functions

- **`EstimatePseudoCell`**  
  Clusters cells and iteratively removes batch effects to estimate pseudo-cells. Useful when actual cell types are unknown.

- **`adj_ortho`**, **`lm.fit.large`**, **`lm.fit.large.new`**  
  Internal helper functions for the orthogonal projection step.

---

## Example Workflow

1. **Load data**  
   ```r
   # Suppose you have:
   #   data.mat: raw counts (cells x genes)
   #   meta.mat: a data.frame with columns 'batch' and 'CellType'
   ```

2. **Check unique batches**  
   ```r
   unique(meta.mat$batch)
   ```

3. **Batch correct**  
   ```r
   result <- BatchCorrect(
     data.mat   = data.mat,
     meta.mat   = meta.mat,
     batch.var  = "batch",
     cell.var   = "CellType"
   )
   # result: a matrix or list(adj, pseudocell=...)
   ```

4. **Downstream analysis**  
   ```r
   # e.g., run PCA, clustering, differential expression on 'result'
   ```

---

## Dependencies

- **R version >= 4.1.1**  
- **Imports**:  
  - [ClusterR](https://cran.r-project.org/package=ClusterR)  
  - [Matrix](https://cran.r-project.org/package=Matrix)  
  - [MASS](https://cran.r-project.org/package=MASS)  
  - [Rfast](https://cran.r-project.org/package=Rfast)  
  - [stats] (base R)  

Install any missing dependencies via `install.packages("...")`.

---

## Contributing

1. **Fork** the repo on GitHub  
2. **Clone** your fork locally  
3. **Create a feature branch** (`git checkout -b new-feature`)  
4. Commit and push your changes  
5. Submit a **Pull Request**

All contributions, bug reports, and suggestions are welcome.

---

## License

This project is licensed under the terms of the **MIT License**.  
See the [LICENSE](./LICENSE) file for details.

---

## Contact
- **Author**: Sujin Seo  
- **Email**: <sujin91kr@gmail.com>  
- If you have any questions or issues, please open an [Issue](https://github.com/sujin91kr/COBRA/issues) on GitHub.

