#' @importFrom stats as.formula lm.fit model.matrix
#' @importFrom Matrix Diagonal sparse.model.matrix
#' @importFrom ClusterR KMeans_arma Optimal_Clusters_KMeans predict_KMeans
#' @keywords internal
"_PACKAGE"

#============================================================
# COBRA package: Batch correction pipeline
#============================================================

#' BatchCorrect
#'
#' Adjusts for batch effects in single-cell RNA-seq data using an orthogonal
#' projection approach. Can also estimate pseudo-cells if cell types are not provided.
#'
#' @param data.mat A matrix of expression data, typically `cells x genes`.
#' @param meta.mat A data.frame containing at least one column for batch (and optionally cell type).
#' @param batch.var A string specifying the batch column in \code{meta.mat}.
#' @param cell.var A string specifying the cell type column in \code{meta.mat}, or \code{NULL} if unknown.
#' @param coi.var A vector of additional covariate column names in \code{meta.mat}.
#' @param stratify Logical. Whether to stratify batch effect correction by cell type.
#'   Default \code{FALSE}.
#' @param reduce Logical. If \code{TRUE}, modifies the input orientation for memory usage.
#'   Default \code{FALSE}.
#' @param sparse.tol Numeric tolerance for zeroing near-zero values, if used. Default \code{1e-6}.
#'
#' @return
#'   - If \code{cell.var} is provided, returns an adjusted expression matrix.
#'   - If \code{cell.var} is \code{NULL}, returns a list with \code{list(adj, pseudocell=...)}.
#' @export
#'
#' @examples
#' \dontrun{
#'   # Suppose data.mat is your count matrix (cells x genes),
#'   # and meta.mat has columns "batch" and "CellType"
#'   corrected <- BatchCorrect(data.mat, meta.mat, batch.var="batch", cell.var="CellType")
#' }
BatchCorrect <- function(data.mat,
                         meta.mat,
                         batch.var = "batch",
                         cell.var = NULL,
                         coi.var   = NULL,
                         stratify  = FALSE,
                         reduce    = FALSE,
                         sparse.tol= 1e-06)
{
  # Re-orient data if needed
  if (reduce) {
    dat <- as.matrix(data.mat)
  } else {
    dat <- Rfast::transpose(as.matrix(data.mat))
  }

  ###############################
  # Cell type not provided --> pseudo-cell estimation
  ###############################
  if (is.null(cell.var)) {
    pseudocell <- EstimatePseudoCell(
      data_mat = data.mat,
      batch    = factor(meta.mat[, batch.var])
    )
    meta.mat$cell <- pseudocell
    orig.cell.var <- NULL
    cell.var <- "cell"
    message("Number of pseudo-cell clusters: ", length(unique(meta.mat$cell)))
  } else {
    orig.cell.var <- cell.var
  }

  ###############################
  # Remove batch effect
  ###############################
  if (!stratify) {
    meta.mat[, batch.var] <- factor(meta.mat[, batch.var])

    # Design formula 만들기
    if (length(unique(meta.mat[, cell.var])) == 1) {
      if (is.null(coi.var)) {
        form <- paste0("~", batch.var)
      } else {
        form <- paste0("~", batch.var, "*(", paste0(coi.var, collapse = "+"), ")")
      }
    } else {
      meta.mat[, cell.var] <- factor(meta.mat[, cell.var])
      if (is.null(coi.var)) {
        form <- paste0("~", batch.var, "*(", cell.var, ")")
      } else {
        form <- paste0("~", batch.var, "*", cell.var, "*(", paste0(coi.var, collapse = "+"), ")")
      }
    }

    design <- sparse.model.matrix(as.formula(form), data = meta.mat)
    batch_loc <- grep(batch.var, colnames(design))

    # Orthogonal Projection
    adj <- adj_ortho(dat, design, batch_loc, reduce)

  } else {
    # stratify = TRUE
    adj <- matrix(NA, nrow(dat), ncol(dat))

    # Identify isolated cell types
    cell_prop <- table(meta.mat[, cell.var], meta.mat[, batch.var]) / c(table(meta.mat[, cell.var]))
    cell_prop_idx <- apply(cell_prop, 1, function(x) sum(x > 0.99))
    cell_even <- names(which(cell_prop_idx == 0))

    sub_cell_loc <- which(meta.mat[, cell.var] %in% cell_even)
    sub_meta <- meta.mat[sub_cell_loc, ]
    sub_dat  <- Rfast::transpose(dat[sub_cell_loc, ])

    # For common cell type
    adj_wo_str <- BatchCorrect(
      data.mat  = sub_dat,
      meta.mat  = sub_meta,
      batch.var = batch.var,
      cell.var  = cell.var,
      stratify  = FALSE
    )
    adj[sub_cell_loc, ] <- adj_wo_str

    # For isolated cell type
    if (length(sub_cell_loc) != nrow(adj)) {
      if (is.null(coi.var)) {
        form <- paste0("~", batch.var, "+", cell.var)
      } else {
        form <- paste0("~", batch.var, "+", cell.var, "+", paste0(coi.var, collapse = "+"))
      }
      design <- sparse.model.matrix(as.formula(form), data = meta.mat)
      batch_loc <- grep(batch.var, colnames(design))
      tmp_adj <- adj_ortho(dat, design, batch_loc, reduce)
      adj[-sub_cell_loc, ] <- tmp_adj[-sub_cell_loc, ]
    }
  }

  # (옵션) re-transpose or apply sparse.tol
  # if (!reduce) {
  #   adj[abs(adj) < sparse.tol] <- 0
  #   adj <- Rfast::transpose(adj)
  #   ...
  # }

  # Return
  if (is.null(orig.cell.var)) {
    return(list(adj = adj, pseudocell = pseudocell))
  } else {
    return(adj)
  }
}


#' adj_ortho: internal helper for orthogonal projection
#'
#' @param dat Transposed or original data (depending on 'reduce')
#' @param design Model matrix (batch and other covariates)
#' @param loc Indices of batch columns in design
#' @param reduce Logical; if TRUE, use different fitting approach
#'
#' @return A batch-effect-adjusted matrix
#' @keywords internal
adj_ortho <- function(dat, design, loc, reduce) {
  n  <- nrow(design)
  xp <- design[, loc]
  xnp <- as.matrix(design[, -loc])

  # Remove any columns in xnp that have only one unique value (besides Intercept)
  uniq_loc <- which(apply(xnp, 2, function(x) length(unique(x))) == 1)
  uniq_loc <- uniq_loc[!grepl('*Intercept*', names(uniq_loc))]
  if (length(uniq_loc) > 0) xnp <- xnp[, -uniq_loc]

  xnpt <- t(xnp)
  xnptxnp <- xnpt %*% xnp
  xnptxnpi <- chol2inv(chol(xnptxnp))

  # xp' (orthogonal projection of xp w.r.t. xnp)
  xpp <- Diagonal(n) %*% xp - xnp %*% (xnptxnpi %*% xnpt %*% xp)
  xpp <- as.matrix(xpp)

  # Fit: choose method by 'reduce'
  if (reduce) {
    # stats::lm.fit() 활용
    coef <- lm.fit(xpp, dat)$coefficients
  } else {
    coef <- lm.fit.large.new(xpp, dat)
  }

  # Remove NA rows if any
  na_loc <- which(is.na(apply(coef, 1, sum)))
  if (length(na_loc) == 0) {
    dat_adj <- dat - xpp %*% coef
  } else {
    dat_adj <- dat - xpp[, -na_loc, drop = FALSE] %*% coef[-na_loc, , drop = FALSE]
  }

  return(dat_adj)
}


#' Large matrix fit (method 1)
#' @keywords internal
lm.fit.large <- function(xpp, dat) {
  xppt <- t(xpp)
  xpptxpp <- as.matrix(xppt %*% xpp)

  if (is.infinite(determinant(xpptxpp)$modulus)) {
    A <- MASS::ginv(xpptxpp) %*% xppt
  } else {
    A <- tryCatch(
      chol2inv(chol(xpptxpp)) %*% xppt,
      error = function(e) MASS::ginv(xpptxpp) %*% xppt
    )
  }

  coef <- matrix(NA, nrow = ncol(xpp), ncol = ncol(dat))
  for (i in seq_len(ncol(dat))) {
    coef[, i] <- A %*% dat[, i]
  }
  return(coef)
}


#' Large matrix fit (method 2, with Rfast transpose)
#' @keywords internal
lm.fit.large.new <- function(xpp, dat) {
  xppt <- Rfast::transpose(xpp)
  xpptxpp <- xppt %*% xpp

  if (is.infinite(determinant(xpptxpp)$modulus)) {
    A <- MASS::ginv(xpptxpp) %*% xppt
  } else {
    A <- tryCatch(
      chol2inv(chol(xpptxpp)) %*% xppt,
      error = function(e) MASS::ginv(xpptxpp) %*% xppt
    )
  }

  coef <- matrix(NA, nrow = ncol(xpp), ncol = ncol(dat))
  for (i in seq_len(ncol(dat))) {
    coef[, i] <- A %*% dat[, i]
  }
  return(coef)
}


#' Estimate Pseudo Cells
#'
#' When \code{cell.var} is NULL, we cluster the data to estimate pseudo-cells
#' and iteratively remove batch effects.
#'
#' @param data_mat Expression matrix (cells x features)
#' @param batch A factor or vector indicating batch labels
#' @param ncl_min Minimum number of clusters to test
#' @param ncl_max Maximum number of clusters to test
#' @param seed Random seed
#' @param r_random_cell Proportion of cells for random sampling
#' @param n_random_cell Absolute minimum number of cells for sampling
#' @param pheno Additional covariates (if any)
#'
#' @return A factor of cluster memberships (pseudo-cells)
#' @export
EstimatePseudoCell <- function(data_mat,
                               batch,
                               ncl_min = 2,
                               ncl_max = 10,
                               seed = 1,
                               r_random_cell = 0.01,
                               n_random_cell = 1000,
                               pheno = NULL)
{
  set.seed(seed)

  # Randomly select some cells
  cell_loc <- sample(seq_len(nrow(data_mat)))[1:max(nrow(data_mat)*r_random_cell, n_random_cell)]
  dat_pcs <- data_mat
  cat("Estimate Pseudo-cells\n")

  # Initial cluster range search
  opt_km <- Optimal_Clusters_KMeans(dat_pcs[cell_loc, ],
                                    max_clusters = c(ncl_min:ncl_max),
                                    plot_clusters = FALSE,
                                    verbose       = FALSE,
                                    criterion     = "silhouette")
  ncl <- max(which(opt_km == max(opt_km))) + ncl_min - 1

  # K-means on full data
  km <- KMeans_arma(dat_pcs, ncl, n_iter = 100, "random_subset", verbose = FALSE)
  cluster <- as.factor(predict_KMeans(dat_pcs, km))

  # BIC for initial
  bic_km <- Optimal_Clusters_KMeans(dat_pcs,
                                    max_clusters = c(ncl_min:ncl_max),
                                    plot_clusters = FALSE,
                                    verbose       = FALSE,
                                    criterion     = "BIC")
  iter <- 1

  while (iter < 20) {
    # Check batch balance
    batch_prop <- table(cluster, batch) / c(table(cluster))
    check_uneven <- ifelse(batch_prop > 0.9, 1, 0)

    # Decide design formula
    if (sum(check_uneven == 1) > length(check_uneven)/length(unique(batch)) * 0.5) {
      if (is.null(unlist(pheno))) {
        design <- model.matrix(~ batch)
      } else {
        design <- model.matrix(~ batch + pheno)
      }
    } else {
      if (is.null(unlist(pheno))) {
        design <- model.matrix(~ batch + cluster)
      } else {
        design <- model.matrix(~ batch + cluster + pheno)
      }
    }
    batch_loc <- grep("^batch", colnames(design))

    # Adjust data
    if (iter == 1) {
      adj_pcs <- adj_ortho(dat_pcs, design, batch_loc, reduce = TRUE)
    } else {
      adj_pcs <- adj_ortho(adj_pcs, design, batch_loc, reduce = TRUE)
    }
    adj_pcs <- as.matrix(adj_pcs)

    # Re-run K-means on adjusted data
    opt_km_adj <- Optimal_Clusters_KMeans(adj_pcs[cell_loc, ],
                                          max_clusters = c(ncl_min:ncl_max),
                                          plot_clusters = FALSE,
                                          verbose       = FALSE,
                                          criterion     = "silhouette")
    ncl_adj <- max(which(opt_km_adj == max(opt_km_adj))) + ncl_min - 1

    km_adj <- KMeans_arma(adj_pcs, ncl_adj, n_iter = 100, "random_subset", verbose = FALSE)
    cluster_adj <- as.factor(predict_KMeans(adj_pcs, km_adj))

    bic_km_adj <- Optimal_Clusters_KMeans(adj_pcs,
                                          max_clusters = c(ncl_min:ncl_max),
                                          plot_clusters = FALSE,
                                          verbose       = FALSE,
                                          criterion     = "BIC")

    if (bic_km_adj[ncl_adj - ncl_min + 1] >= bic_km[ncl_adj - ncl_min + 1]) {
      break
    } else {
      message("Iteration ", iter, ": ", ncl_adj, " clusters")
      iter <- iter + 1
      bic_km <- bic_km_adj
      cluster <- cluster_adj
    }
  }

  return(cluster)
}
