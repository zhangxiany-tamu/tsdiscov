#' Covariance matrix features for multivariate time series
#'
#' Extracts features based on the covariance matrix of the time series.
#' Unlike correlation features which use standardized data, these use the
#' actual variances and covariances.
#'
#' @param X Matrix (N x T) where N is number of series, T is timepoints
#' @return Named list of 7 covariance-based features
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

  # Feature 2: Determinant
  det_Sigma <- tryCatch(
    det(Sigma),
    error = function(e) NA
  )
  cov_determinant <- det_Sigma

  # Feature 3: Log determinant (generalized variance)
  cov_log_determinant <- if(!is.na(det_Sigma) && det_Sigma > 1e-10) {
    log(det_Sigma)
  } else {
    NA
  }

  # Feature 4: Condition number
  lambda_min <- eigenvalues[N]
  lambda_max <- eigenvalues[1]
  cov_condition_number <- if(lambda_min > 1e-10) {
    lambda_max / lambda_min
  } else {
    NA
  }

  # Feature 5: Frobenius norm
  # ||Σ||_F = sqrt(sum(Σ^2))
  cov_frobenius_norm <- sqrt(sum(Sigma^2))

  # Feature 6: Spectral norm (largest eigenvalue)
  cov_spectral_norm <- lambda_max

  # Return named list (removed cov_nuclear_norm - always equals cov_trace)
  list(
    cov_trace = cov_trace,
    cov_determinant = cov_determinant,
    cov_log_determinant = cov_log_determinant,
    cov_condition_number = cov_condition_number,
    cov_frobenius_norm = cov_frobenius_norm,
    cov_spectral_norm = cov_spectral_norm
  )
}
