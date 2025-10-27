#' Covariance matrix features for multivariate time series
#'
#' Extracts features based on the covariance matrix of the time series.
#' Unlike correlation features which use standardized data, these use the
#' actual variances and covariances.
#'
#' @param X Matrix (N x T) where N is number of series, T is timepoints
#' @return Named list of 5 covariance-based features
#' @export
#' @examples
#' \dontrun{
#' # 3 time series with different variances
#' X <- matrix(0, nrow = 3, ncol = 150)
#' X[1,] <- rnorm(150, sd = 1)
#' X[2,] <- rnorm(150, sd = 2)
#' X[3,] <- rnorm(150, sd = 0.5)
#' cov_features <- ts_mv_covariance(X)
#' }
ts_mv_covariance <- function(X) {
  N <- nrow(X)
  T_len <- ncol(X)

  # Compute covariance matrix
  Sigma <- cov(t(X))

  # Eigenvalue decomposition
  eigen_result <- eigen(Sigma, symmetric = TRUE)
  eigenvalues <- eigen_result$values

  # Ensure non-negative (numerical precision)
  eigenvalues <- pmax(eigenvalues, 0)

  # Feature 1: Trace (total variance across all series)
  cov_trace <- sum(diag(Sigma))

  # Feature 2: Log determinant (robust, using eigenvalue floor)
  cov_log_determinant <- robust_log_det(eigenvalues, floor = 1e-10)

  # Feature 3: Condition number (robust, returns NA for near-singular)
  cov_condition_number <- robust_cond_number(eigenvalues, eps = 1e-10)

  # Feature 4: Frobenius norm
  # ||Σ||_F = sqrt(sum(Σ^2))
  cov_frobenius_norm <- sqrt(sum(Sigma^2))

  # Feature 5: Spectral norm (largest eigenvalue)
  cov_spectral_norm <- max(eigenvalues)

  # Return named list (now 5 features instead of 6)
  list(
    cov_trace = cov_trace,
    cov_log_determinant = cov_log_determinant,
    cov_condition_number = cov_condition_number,
    cov_frobenius_norm = cov_frobenius_norm,
    cov_spectral_norm = cov_spectral_norm
  )
}
