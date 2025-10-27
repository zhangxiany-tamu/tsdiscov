#' Synchronization features for multivariate time series
#'
#' Extracts features characterizing time-domain synchronization between series.
#' Uses cross-correlation functions (CCF) and mutual information to measure
#' how synchronized the series are.
#'
#' @param X Matrix (N x T) where N is number of series, T is timepoints
#' @param max_lag Integer; maximum lag for CCF computation (default 20).
#'   Will be capped at floor(T/4) where T is the series length.
#' @param bins Integer; number of bins for mutual information computation (default 10)
#' @param ... Additional arguments (currently ignored)
#' @return Named list of 8 synchronization-based features
#' @export
#' @examples
#' \dontrun{
#' # Create synchronized signals
#' t <- seq(0, 10, length.out = 200)
#' X <- matrix(0, nrow = 3, ncol = 200)
#' X[1,] <- sin(2*pi*t) + rnorm(200, sd=0.1)
#' X[2,] <- sin(2*pi*t + 0.2) + rnorm(200, sd=0.1)  # Small phase lag
#' X[3,] <- sin(2*pi*t) + rnorm(200, sd=0.1)
#' sync_features <- ts_mv_synchronization(X)
#' }
ts_mv_synchronization <- function(X, max_lag = 20, bins = 10, ...) {
  N <- nrow(X)
  T_len <- ncol(X)

  if (N < 2) {
    # Single series - no synchronization
    return(list(
      ccf_max_mean = NA, ccf_max_median = NA, ccf_max_std = NA,
      ccf_lag_mean = NA, ccf_lag_std = NA,
      mutual_info_mean = NA, mutual_info_max = NA, mutual_info_std = NA
    ))
  }

  # Cap max_lag at floor(T/4) to ensure reasonable CCF computation
  max_lag <- min(max_lag, floor(T_len / 4))

  # Compute CCF for all pairs using Rcpp (much faster)
  ccf_result <- cpp_pairwise_ccf(X, max_lag)
  ccf_max_vals <- ccf_result$ccf_max_vals
  ccf_lag_vals <- ccf_result$ccf_lag_vals
  ccf_zero_vals <- ccf_result$ccf_zero_vals

  # Compute mutual information for all pairs using Rcpp
  mi_vals <- cpp_pairwise_mi(X, bins = bins)

  # Feature 1: Mean of maximum CCF
  ccf_max_mean <- mean(ccf_max_vals, na.rm = TRUE)

  # Feature 2: Median of maximum CCF
  ccf_max_median <- median(ccf_max_vals, na.rm = TRUE)

  # Feature 3: Std of maximum CCF
  ccf_max_std <- sd(ccf_max_vals, na.rm = TRUE)

  # Feature 4: Mean lag at maximum CCF (removed ccf_sync_frac - redundant with corr_strong_frac)
  ccf_lag_mean <- mean(abs(ccf_lag_vals), na.rm = TRUE)

  # Feature 5: Std of lags
  ccf_lag_std <- sd(ccf_lag_vals, na.rm = TRUE)

  # Feature 6: Mean mutual information (removed ccf_zero_lag_mean - redundant with corr_mean)
  mutual_info_mean <- mean(mi_vals, na.rm = TRUE)

  # Feature 7: Maximum mutual information
  mutual_info_max <- max(mi_vals, na.rm = TRUE)

  # Feature 8: Std of mutual information
  mutual_info_std <- sd(mi_vals, na.rm = TRUE)

  # Return named list (removed ccf_sync_frac and ccf_zero_lag_mean - redundant)
  list(
    ccf_max_mean = ccf_max_mean,
    ccf_max_median = ccf_max_median,
    ccf_max_std = ccf_max_std,
    ccf_lag_mean = ccf_lag_mean,
    ccf_lag_std = ccf_lag_std,
    mutual_info_mean = mutual_info_mean,
    mutual_info_max = mutual_info_max,
    mutual_info_std = mutual_info_std
  )
}

#' Compute mutual information between two vectors
#'
#' Simple histogram-based estimator of mutual information
#'
#' @param x Numeric vector
#' @param y Numeric vector
#' @param bins Number of bins for histogram (default 10)
#' @return Mutual information estimate
#' @keywords internal
compute_mutual_information <- function(x, y, bins = 10) {
  # Remove NA
  valid <- !is.na(x) & !is.na(y)
  x <- x[valid]
  y <- y[valid]

  if (length(x) < 10) return(NA)

  # Discretize into bins
  x_breaks <- seq(min(x), max(x), length.out = bins + 1)
  y_breaks <- seq(min(y), max(y), length.out = bins + 1)

  x_binned <- cut(x, breaks = x_breaks, include.lowest = TRUE, labels = FALSE)
  y_binned <- cut(y, breaks = y_breaks, include.lowest = TRUE, labels = FALSE)

  # Joint histogram
  joint_hist <- table(x_binned, y_binned)
  p_xy <- joint_hist / sum(joint_hist)

  # Marginal histograms
  p_x <- rowSums(p_xy)
  p_y <- colSums(p_xy)

  # Mutual information: I(X;Y) = sum p(x,y) * log(p(x,y) / (p(x)*p(y)))
  mi <- 0
  for (i in 1:bins) {
    for (j in 1:bins) {
      if (p_xy[i,j] > 0 && p_x[i] > 0 && p_y[j] > 0) {
        mi <- mi + p_xy[i,j] * log(p_xy[i,j] / (p_x[i] * p_y[j]))
      }
    }
  }

  return(mi)
}
