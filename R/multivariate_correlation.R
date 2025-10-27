#' Correlation structure features for multivariate time series
#'
#' Extracts features characterizing the correlation structure between time series.
#'
#' @param X Matrix (N x T) where N is number of series, T is timepoints
#' @param method Character; correlation method. One of "pearson" (default),
#'   "spearman", or "kendall"
#' @param ... Additional arguments (currently ignored)
#' @return Named list of 15 correlation-based features
#' @export
#' @examples
#' \dontrun{
#' # 4 time series with varying correlations
#' X <- matrix(rnorm(4 * 150), nrow = 4)
#' X[2,] <- 0.8 * X[1,] + 0.2 * rnorm(150)
#' corr_features <- ts_mv_correlation(X)
#' }
ts_mv_correlation <- function(X, method = "pearson", ...) {
  # Validate method
  method <- match.arg(method, c("pearson", "spearman", "kendall"))

  N <- nrow(X)
  T_len <- ncol(X)

  # Compute correlation matrix with specified method
  R <- cor(t(X), method = method)

  # Extract off-diagonal elements (the actual correlations between series)
  if (N >= 2) {
    off_diag <- R[lower.tri(R)]
    n_pairs <- length(off_diag)
  } else {
    # Single series - no correlations
    return(list(
      corr_mean = NA, corr_median = NA, corr_std = NA,
      corr_max = NA, corr_min = NA, corr_range = NA, corr_iqr = NA,
      corr_skewness = NA, corr_positive_frac = NA, corr_strong_frac = NA,
      corr_weak_frac = NA, corr_spectral_radius = NA, corr_determinant = NA,
      corr_log_determinant = NA, corr_frobenius_norm = NA
    ))
  }

  # Feature 1: Mean correlation
  corr_mean <- mean(off_diag)

  # Feature 2: Median correlation
  corr_median <- median(off_diag)

  # Feature 3: Standard deviation of correlations
  corr_std <- if (length(off_diag) > 1) sd(off_diag) else 0

  # Feature 4: Maximum correlation
  corr_max <- max(off_diag)

  # Feature 5: Minimum correlation
  corr_min <- min(off_diag)

  # Feature 6: Range of correlations
  corr_range <- corr_max - corr_min

  # Feature 7: IQR of correlations
  corr_iqr <- IQR(off_diag)

  # Feature 8: Skewness of correlation distribution
  if (n_pairs >= 3) {
    m <- mean(off_diag)
    s <- sd(off_diag)
    if (s > 1e-10) {
      corr_skewness <- mean((off_diag - m)^3) / s^3
    } else {
      corr_skewness <- 0
    }
  } else {
    corr_skewness <- NA
  }

  # Feature 9: Fraction of positive correlations
  corr_positive_frac <- mean(off_diag > 0)

  # Feature 10: Fraction of strong correlations (|r| > 0.5)
  corr_strong_frac <- mean(abs(off_diag) > 0.5)

  # Feature 11: Fraction of weak correlations (|r| < 0.1)
  corr_weak_frac <- mean(abs(off_diag) < 0.1)

  # Feature 12: Spectral radius (largest eigenvalue)
  eigen_vals <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  corr_spectral_radius <- max(eigen_vals)

  # Feature 13: Log determinant (robust, using eigenvalue floor)
  corr_log_determinant <- robust_log_det(eigen_vals, floor = 1e-10)

  # Feature 14: Condition number (robust, returns NA for near-singular)
  corr_condition_number <- robust_cond_number(eigen_vals, eps = 1e-10)

  # Feature 15: Frobenius norm of (R - I)
  # ||R - I||_F = sqrt(sum((R - I)^2))
  R_centered <- R - diag(N)
  corr_frobenius_norm <- sqrt(sum(R_centered^2))

  # Return named list
  list(
    corr_mean = corr_mean,
    corr_median = corr_median,
    corr_std = corr_std,
    corr_max = corr_max,
    corr_min = corr_min,
    corr_range = corr_range,
    corr_iqr = corr_iqr,
    corr_skewness = corr_skewness,
    corr_positive_frac = corr_positive_frac,
    corr_strong_frac = corr_strong_frac,
    corr_weak_frac = corr_weak_frac,
    corr_spectral_radius = corr_spectral_radius,
    corr_log_determinant = corr_log_determinant,
    corr_condition_number = corr_condition_number,
    corr_frobenius_norm = corr_frobenius_norm
  )
}
