# This file contains helper functions for BayesCNet

#' Find overlapping genomic coordinates
#'
#' @param peak_list Character vector of peaks in format "chr_start_end"
#' @param target_region Target region in same format
#'
#' @return Character vector of overlapping peaks
#' @keywords internal
find_overlapping_coordinates <- function(peak_list, target_region) {
  # Parse target region
  target_parts <- strsplit(target_region, "_")[[1]]
  target_chr <- target_parts[1]
  target_start <- as.numeric(target_parts[2])
  target_end <- as.numeric(target_parts[3])

  # Find overlapping peaks
  overlapping <- sapply(peak_list, function(peak) {
    parts <- strsplit(peak, "_")[[1]]
    chr <- parts[1]
    start <- as.numeric(parts[2])
    end <- as.numeric(parts[3])

    # Check if same chromosome and overlapping
    if (chr == target_chr &&
        ((start >= target_start && start <= target_end) ||
         (end >= target_start && end <= target_end) ||
         (start <= target_start && end >= target_end))) {
      return(TRUE)
    }
    return(FALSE)
  })

  return(peak_list[overlapping])
}

#' Estimate size factors for matrix normalization
#'
#' @param mat Count matrix
#'
#' @return Numeric vector of size factors
#' @keywords internal
estimateSizeFactorsForMatrix <- function(mat) {
  # Simple library size normalization
  # Can be replaced with more sophisticated methods
  col_sums <- colSums(mat)
  size_factors <- col_sums / median(col_sums)
  return(size_factors)
}

#' Process a single gene for Bayesian inference
#'
#' @param gene_idx Gene index
#' @param gene_info Data frame with gene information
#' @param data_rna RNA data matrix
#' @param data_atac ATAC data matrix
#' @param cell_types_vector Cell type assignments
#' @param Sigma2_est Estimated covariance matrix
#' @param window Regulatory window size
#' @param regularization Regularization parameter
#'
#' @return Data frame with inference results
#' @keywords internal
process_gene_inference <- function(
    gene_idx,
    gene_info,
    data_rna,
    data_atac,
    cell_types_vector,
    Sigma2_est,
    window = 250000,
    regularization = 1e-6
) {

  # Extract gene information
  gene_name <- gene_info$gene[gene_idx]
  chr <- gene_info$chr[gene_idx]
  tss <- gene_info$tss[gene_idx]

  # Define promoter and regulatory regions
  promoter <- paste0(chr, "_", max(tss - 500, 1), "_", tss)
  reg_start <- max(tss - window, 1)
  reg_end <- tss + window
  regulatory_region <- paste0(chr, "_", reg_start, "_", reg_end)

  # Get peaks in regulatory region
  potential_peaks <- rownames(data_atac)
  promoter_peaks <- find_overlapping_coordinates(potential_peaks, promoter)
  if (length(promoter_peaks) != 1) {
    promoter_peaks <- promoter
  }

  enhancer_peaks <- find_overlapping_coordinates(potential_peaks, regulatory_region)
  enhancer_peaks <- setdiff(enhancer_peaks, promoter_peaks)

  if (length(enhancer_peaks) < 2) {
    return(NULL)
  }

  # Build data matrices
  # Check if gene exists in RNA data
  if (!gene_name %in% rownames(data_rna)) {
    return(NULL)
  }

  gene_expr <- data_rna[gene_name, ]
  peak_idx <- match(enhancer_peaks, potential_peaks)
  peak_idx <- peak_idx[!is.na(peak_idx)]

  X <- data_atac[peak_idx, , drop = FALSE]
  X <- as.matrix(X)
  tX <- t(X)

  # Split by cell type
  X_list <- lapply(split(seq_len(nrow(tX)),
                         factor(cell_types_vector, levels = unique(cell_types_vector))),
                   function(idx) tX[idx, , drop = FALSE])

  # Create block matrix
  X_block <- Matrix::bdiag(X_list)
  X_block <- as.matrix(X_block)

  # Prepare Y
  Y <- as.matrix(gene_expr)

  # Get dimensions
  celltypes_order <- unique(cell_types_vector)
  C <- length(celltypes_order)
  P <- ncol(tX)

  # Reorder Sigma2_est to match cell type order
  Sigma2_est <- Sigma2_est[celltypes_order, celltypes_order]

  # Compute inference using hierarchical Bayesian model
  Sigma2_est_sqrt <- expm::sqrtm(Sigma2_est)
  SigIp <- kronecker(Sigma2_est_sqrt, diag(P))
  tilde_X <- X_block %*% SigIp

  Z <- do.call(rbind, X_list)

  # OLS for initialization
  tryCatch({
    ols_fit <- lm(Y ~ Z - 1)
    Y_hat_OLS <- predict(ols_fit)
    residuals <- Y - Y_hat_OLS
    sigma_hat_sq <- var(residuals)

    # SURE estimator for adaptive shrinkage
    numerator <- sum((Y_hat_OLS - Y)^2)
    denominator <- P * sigma_hat_sq
    w_SURE <- as.numeric((numerator / denominator) - 1)

    # Ensure w_SURE is positive
    w_SURE <- max(w_SURE, 0.01)

  }, error = function(e) {
    # If OLS fails, use default values
    sigma_hat_sq <- var(Y)
    w_SURE <- 1
  })

  # Prior covariance
  Sigma_1_est <- w_SURE * solve(t(Z) %*% Z + regularization * diag(ncol(Z)))
  V_0 <- kronecker(diag(C), Sigma_1_est)

  # Prior mean (zero)
  g_0_rep <- rep(0, P * C)

  # Posterior computation
  V_n <- solve(solve(V_0) + t(tilde_X) %*% tilde_X)
  g_n <- V_n %*% (solve(V_0) %*% g_0_rep + t(tilde_X) %*% Y)

  # Hyperparameters for inverse gamma prior on noise variance
  alpha_0 <- 1
  delta_0 <- 1
  n_total <- length(Y)
  alpha_n <- alpha_0 + n_total / 2
  delta_n <- delta_0 + 0.5 * (t(Y) %*% Y + t(g_0_rep) %*% solve(V_0) %*% g_0_rep -
                                t(g_n) %*% solve(V_n) %*% g_n)
  mean_sigma2 <- as.numeric(delta_n / (alpha_n - 1))

  # Extract posterior means and variances
  B_posterior_mean <- SigIp %*% g_n
  B_posterior_mean <- matrix(B_posterior_mean, nrow = P, ncol = C, byrow = FALSE)

  B_posterior_covar <- SigIp %*% (mean_sigma2 * V_n) %*% t(SigIp)
  B_posterior_var <- diag(B_posterior_covar)
  B_posterior_var <- matrix(B_posterior_var, nrow = P, ncol = C, byrow = FALSE)

  # Format results
  results_all <- NULL
  for (j in 1:C) {
    conns <- data.frame(
      CellType = celltypes_order[j],
      Gene = gene_name,
      Peak1 = promoter_peaks[1],
      Peak2 = enhancer_peaks,
      Pmean = B_posterior_mean[, j],
      Pvar = B_posterior_var[, j],
      Importance = abs(B_posterior_mean[, j]) / sqrt(B_posterior_var[, j]),
      stringsAsFactors = FALSE
    )
    results_all <- rbind(results_all, conns)
  }

  return(results_all)
}

#' Aggregate a single cell type
#'
#' @param rna_data RNA count matrix for one cell type
#' @param atac_data ATAC count matrix for one cell type
#' @param cell_coord Cell coordinates for this cell type
#' @param k_neigh Number of neighbors
#' @param max_overlap Maximum overlap allowed
#' @param seed Random seed
#' @param verbose Print messages
#'
#' @return List with aggregated data
#' @keywords internal
aggregate_single_celltype <- function(rna_data, atac_data, cell_coord,
                                      k_neigh = 50, max_overlap = 0.8,
                                      seed = 123, verbose = FALSE) {

  n_cells <- nrow(cell_coord)

  if (n_cells <= k_neigh) {
    # If fewer cells than k, aggregate all cells into one metacell
    if (verbose) message("    Fewer cells than k, creating single metacell")

    rna_agg <- rowSums(rna_data)
    rna_agg <- matrix(rna_agg, ncol = 1)
    rownames(rna_agg) <- rownames(rna_data)

    atac_agg <- rowSums(atac_data > 0)  # Binarize ATAC
    atac_agg <- matrix(atac_agg, ncol = 1)
    rownames(atac_agg) <- rownames(atac_data)

    cell_sample <- matrix(1:n_cells, nrow = 1)

  } else {
    # KNN aggregation
    if (verbose) message("    Running KNN aggregation")

    # Find k nearest neighbors
    nn_map <- as.data.frame(FNN::knn.index(cell_coord, k = (k_neigh - 1)))
    rownames(nn_map) <- rownames(cell_coord)
    nn_map$agg_cell <- 1:nrow(nn_map)
    good_choices <- 1:nrow(nn_map)

    # Sample cells with minimal overlap
    set.seed(seed)
    choice <- sample(good_choices, size = 1)
    chosen <- good_choices[choice]
    good_choices <- good_choices[-choice]

    it <- 0
    max_it <- n_cells / ((1 - max_overlap) * k_neigh)

    while (length(good_choices) > 0 && it < max_it) {
      it <- it + 1
      choice <- sample(seq_along(good_choices), size = 1)
      new_chosen <- c(chosen, good_choices[choice])
      good_choices <- good_choices[-choice]

      # Get current cell samples
      cell_sample_temp <- nn_map[new_chosen, -ncol(nn_map)]

      # Check overlap
      if (nrow(cell_sample_temp) > 1) {
        combs <- expand.grid(1:(length(new_chosen)-1), length(new_chosen))
        shared <- apply(combs, 1, function(x) {
          cells1 <- unlist(cell_sample_temp[x[1], ])
          cells2 <- unlist(cell_sample_temp[x[2], ])
          length(intersect(cells1, cells2))
        })

        if (max(shared) < max_overlap * k_neigh) {
          chosen <- new_chosen
        }
      } else {
        chosen <- new_chosen
      }
    }

    # Aggregate selected cells
    cell_sample <- as.matrix(nn_map[chosen, -ncol(nn_map)])
    n_metacells <- length(chosen)

    # Aggregate RNA
    rna_agg <- matrix(0, nrow = nrow(rna_data), ncol = n_metacells)
    rownames(rna_agg) <- rownames(rna_data)
    for (i in 1:n_metacells) {
      cell_indices <- c(chosen[i], unlist(cell_sample[i, ]))
      rna_agg[, i] <- rowSums(rna_data[, cell_indices, drop = FALSE])
    }

    # Aggregate ATAC (binarized)
    atac_agg <- matrix(0, nrow = nrow(atac_data), ncol = n_metacells)
    rownames(atac_agg) <- gsub("-", "_", rownames(atac_data))
    for (i in 1:n_metacells) {
      cell_indices <- c(chosen[i], unlist(cell_sample[i, ]))
      atac_agg[, i] <- rowSums(atac_data[, cell_indices, drop = FALSE] > 0)
    }

    # Add the chosen cell itself to the sample matrix
    cell_sample <- cbind(chosen, cell_sample)
  }

  return(list(
    rna = rna_agg,
    atac = atac_agg,
    cell_sample = cell_sample
  ))
}
