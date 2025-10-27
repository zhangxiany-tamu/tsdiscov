#' Lag Structure Features
#'
#' Analyzes the temporal relationships and time delays between series by examining
#' the distribution of optimal lags in cross-correlation functions. These features
#' are useful for detecting leader-follower dynamics, feedback loops, and temporal
#' dependencies in multivariate systems.
#'
#' @param X Numeric matrix (N x T) where rows are time series
#' @param max_lag Integer, maximum lag to consider (default: 20)
#' @param ... Additional arguments (for compatibility)
#'
#' @return Named list with 4 features:
#' \itemize{
#'   \item \code{lag_structure_entropy}: Shannon entropy of lag distribution
#'   \item \code{lag_structure_peak}: Most common optimal lag
#'   \item \code{lag_structure_mean}: Mean of optimal lags (weighted by CCF magnitude)
#'   \item \code{lag_structure_std}: Standard deviation of optimal lags
#' }
#'
#' @details
#' For each pair of series, computes the cross-correlation function (CCF) and
#' identifies the lag that maximizes |CCF|. The distribution of these optimal
#' lags across all pairs provides insight into the temporal structure:
#'
#' \itemize{
#'   \item Low entropy: Consistent temporal relationships (e.g., synchronized or
#'         consistent leader-follower)
#'   \item High entropy: Diverse temporal relationships
#'   \item lag_structure_peak: Dominant temporal relationship
#'   \item Non-zero mean: Asymmetric leader-follower dynamics
#' }
#'
#' @references
#' Bressler, S. L., & Seth, A. K. (2011). Wiener–Granger causality: a well
#' established methodology. Neuroimage, 58(2), 323-329.
#'
#' @export
ts_mv_lag_structure <- function(X, max_lag = 20, ...) {
  N <- nrow(X)
  T_len <- ncol(X)

  # Check for sufficient data
  if (N < 2) {
    return(list(
      lag_structure_entropy = NA_real_,
      lag_structure_peak = NA_real_,
      lag_structure_mean = NA_real_,
      lag_structure_std = NA_real_
    ))
  }

  if (T_len < 10) {
    return(list(
      lag_structure_entropy = NA_real_,
      lag_structure_peak = NA_real_,
      lag_structure_mean = NA_real_,
      lag_structure_std = NA_real_
    ))
  }

  # Cap max_lag at T/4
  max_lag <- min(max_lag, floor(T_len / 4))

  # Compute optimal lags for all pairs
  n_pairs <- N * (N - 1) / 2
  optimal_lags <- numeric(n_pairs)
  ccf_magnitudes <- numeric(n_pairs)

  idx <- 1
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      # Compute CCF
      ccf_result <- ccf(X[i, ], X[j, ], lag.max = max_lag, plot = FALSE)

      # Find lag with maximum absolute CCF
      max_idx <- which.max(abs(ccf_result$acf))
      optimal_lags[idx] <- ccf_result$lag[max_idx]
      ccf_magnitudes[idx] <- abs(ccf_result$acf[max_idx])

      idx <- idx + 1
    }
  }

  # Compute features from lag distribution

  # 1. Entropy of lag distribution (how diverse are the temporal relationships?)
  lag_counts <- table(optimal_lags)
  lag_probs <- lag_counts / sum(lag_counts)
  lag_structure_entropy <- -sum(lag_probs * log(lag_probs + 1e-10))

  # 2. Peak lag (most common temporal relationship)
  lag_structure_peak <- as.numeric(names(lag_counts)[which.max(lag_counts)])

  # 3. Weighted mean lag (weighted by CCF magnitude)
  # Lags with stronger correlations have more weight
  weights <- ccf_magnitudes / sum(ccf_magnitudes + 1e-10)
  lag_structure_mean <- sum(optimal_lags * weights)

  # 4. Standard deviation of lags
  lag_structure_std <- sd(optimal_lags)

  list(
    lag_structure_entropy = lag_structure_entropy,
    lag_structure_peak = lag_structure_peak,
    lag_structure_mean = lag_structure_mean,
    lag_structure_std = lag_structure_std
  )
}
