#' PCA-based multivariate features
#'
#' Extracts features based on principal component analysis of the time series.
#' Performs PCA on the correlation matrix of the N time series.
#'
#' @param X Matrix (N x T) where N is number of series, T is timepoints
#' @return Named list of 15 PCA-based features
#' @export
#' @examples
#' \dontrun{
#' # 5 correlated time series
#' X <- matrix(rnorm(5 * 200), nrow = 5)
#' X[2,] <- 0.7 * X[1,] + 0.3 * rnorm(200)
#' pca_features <- ts_mv_pca(X)
#' }
ts_mv_pca <- function(X) {
  N <- nrow(X)
  T_len <- ncol(X)

  # Compute correlation matrix
  R <- cor(t(X))

  # Eigenvalue decomposition
  eigen_result <- eigen(R, symmetric = TRUE)
  eigenvalues <- eigen_result$values

  # Ensure non-negative (numerical precision)
  eigenvalues <- pmax(eigenvalues, 0)

  # Normalize eigenvalues to sum to N (for correlation matrix, trace = N)
  eigenvalues <- eigenvalues / sum(eigenvalues) * N

  # Variance explained (as proportion)
  var_explained <- eigenvalues / sum(eigenvalues)
  cumvar_explained <- cumsum(var_explained)

  # Feature 1-3: Variance explained by first 3 PCs
  pca_var_pc1 <- var_explained[1]
  pca_var_pc2 <- if(N >= 2) var_explained[2] else NA
  pca_var_pc3 <- if(N >= 3) var_explained[3] else NA

  # Feature 4: Cumulative variance of first 3 PCs
  pca_var_cumsum_3 <- if(N >= 3) cumvar_explained[3] else cumvar_explained[min(N, 3)]

  # Feature 5: Ratio of PC1 to PC2 variance
  pca_var_ratio_1_2 <- if(N >= 2 && eigenvalues[2] > 1e-10) {
    eigenvalues[1] / eigenvalues[2]
  } else {
    NA
  }

  # Feature 6: Effective rank (participation ratio)
  # PR = (sum λ_i)^2 / sum(λ_i^2)
  pca_effective_rank <- sum(eigenvalues)^2 / sum(eigenvalues^2)

  # Feature 7: Participation ratio (inverse)
  pca_participation_ratio <- sum(eigenvalues^2) / sum(eigenvalues)^2

  # Feature 8: Eigenvalue entropy
  # H = -sum(p_i * log(p_i)) where p_i = λ_i / sum(λ_i)
  p <- eigenvalues / sum(eigenvalues)
  p <- p[p > 1e-10]  # Avoid log(0)
  pca_entropy <- -sum(p * log(p))

  # Feature 9: Stable rank
  # SR = sum(λ_i) / λ_max
  pca_stable_rank <- sum(eigenvalues) / eigenvalues[1]

  # Feature 10: Condition number
  lambda_min <- eigenvalues[N]
  pca_condition_number <- if(lambda_min > 1e-10) {
    eigenvalues[1] / lambda_min
  } else {
    NA
  }

  # Feature 11-12: Number of PCs needed for 90% and 95% variance
  pca_num_pcs_90pct <- which(cumvar_explained >= 0.90)[1]
  pca_num_pcs_95pct <- which(cumvar_explained >= 0.95)[1]

  # Feature 13: Maximum eigenvalue gap
  if (N >= 2) {
    gaps <- diff(eigenvalues)
    pca_eigenvalue_gap_max <- max(abs(gaps))
  } else {
    pca_eigenvalue_gap_max <- NA
  }

  # Feature 14: Eigenvalue decay rate
  # Fit log(eigenvalue) ~ rank to get decay slope
  if (N >= 3) {
    positive_eigs <- eigenvalues[eigenvalues > 1e-10]
    if (length(positive_eigs) >= 2) {
      ranks <- 1:length(positive_eigs)
      fit <- tryCatch(
        lm(log(positive_eigs) ~ ranks),
        error = function(e) NULL
      )
      pca_eigenvalue_decay <- if(!is.null(fit)) coef(fit)[2] else NA
    } else {
      pca_eigenvalue_decay <- NA
    }
  } else {
    pca_eigenvalue_decay <- NA
  }

  # Feature 15: Scree elbow detection
  # Use maximum second derivative as elbow point
  if (N >= 4) {
    second_diff <- diff(diff(eigenvalues))
    pca_scree_elbow <- which.max(abs(second_diff)) + 1
  } else {
    pca_scree_elbow <- NA
  }

  # Return named list
  list(
    pca_var_pc1 = pca_var_pc1,
    pca_var_pc2 = pca_var_pc2,
    pca_var_pc3 = pca_var_pc3,
    pca_var_cumsum_3 = pca_var_cumsum_3,
    pca_var_ratio_1_2 = pca_var_ratio_1_2,
    pca_effective_rank = pca_effective_rank,
    pca_participation_ratio = pca_participation_ratio,
    pca_entropy = pca_entropy,
    pca_stable_rank = pca_stable_rank,
    pca_condition_number = pca_condition_number,
    pca_num_pcs_90pct = pca_num_pcs_90pct,
    pca_num_pcs_95pct = pca_num_pcs_95pct,
    pca_eigenvalue_gap_max = pca_eigenvalue_gap_max,
    pca_eigenvalue_decay = pca_eigenvalue_decay,
    pca_scree_elbow = pca_scree_elbow
  )
}
