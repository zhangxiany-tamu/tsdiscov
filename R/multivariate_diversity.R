#' Diversity and complexity features for multivariate time series
#'
#' Extracts features characterizing the diversity and complexity of the
#' multivariate system. Includes entropy-based measures and heterogeneity metrics.
#'
#' @param X Matrix (N x T) where N is number of series, T is timepoints
#' @param bins Integer; number of bins for Shannon entropy computation (default 10)
#' @param ... Additional arguments (currently ignored)
#' @return Named list of 7 diversity/complexity features
#' @export
#' @examples
#' \dontrun{
#' # Heterogeneous system
#' X <- matrix(0, nrow = 4, ncol = 200)
#' X[1,] <- arima.sim(list(ar = 0.9), 200)  # Strong AR
#' X[2,] <- arima.sim(list(ar = 0.3), 200)  # Weak AR
#' X[3,] <- rnorm(200)  # White noise
#' X[4,] <- cumsum(rnorm(200))  # Random walk
#' diversity_features <- ts_mv_diversity(X)
#' }
ts_mv_diversity <- function(X, bins = 10, ...) {
  N <- nrow(X)
  T_len <- ncol(X)

  if (N < 2) {
    # Single series - no diversity
    return(list(
      mean_diversity = NA, acf1_diversity = NA,
      complexity_variance = NA, entropy_diversity = NA,
      range_diversity = NA, skewness_diversity = NA,
      kurtosis_diversity = NA
    ))
  }

  # Compute basic statistics for each series using Rcpp (fast)
  stats <- cpp_series_statistics(X)
  means <- stats$means
  variances <- stats$variances
  ranges <- stats$ranges
  skewness_vals <- stats$skewness
  kurtosis_vals <- stats$kurtosis

  # ACF lag 1 for each series (still using R's acf for now)
  acf1_vals <- numeric(N)
  for (i in 1:N) {
    acf_result <- tryCatch({
      acf(X[i,], lag.max = 1, plot = FALSE, na.action = na.pass)
    }, error = function(e) NULL)

    if (!is.null(acf_result) && length(acf_result$acf) >= 2) {
      acf1_vals[i] <- acf_result$acf[2]
    } else {
      acf1_vals[i] <- NA
    }
  }

  # Shannon entropy for each series using Rcpp (using provided bins parameter)
  entropy_vals <- numeric(N)
  for (i in 1:N) {
    entropy_vals[i] <- tryCatch({
      cpp_mv_shannon_entropy(X[i,], bins = bins)
    }, error = function(e) NA)
  }

  # Feature 1: Mean diversity (coefficient of variation of means)
  # (removed variance_diversity - always 0 when standardized)
  mean_diversity <- if (mean(abs(means), na.rm = TRUE) > 1e-10) {
    sd(means, na.rm = TRUE) / (mean(abs(means), na.rm = TRUE) + 1e-10)
  } else {
    sd(means, na.rm = TRUE)
  }

  # Feature 2: ACF diversity (std of ACF lag-1 values)
  acf1_diversity <- sd(acf1_vals, na.rm = TRUE)

  # Feature 3: Complexity variance (variance of entropies)
  complexity_variance <- var(entropy_vals, na.rm = TRUE)

  # Feature 4: Entropy diversity (coefficient of variation of entropies)
  entropy_diversity <- if (mean(entropy_vals, na.rm = TRUE) > 1e-10) {
    sd(entropy_vals, na.rm = TRUE) / mean(entropy_vals, na.rm = TRUE)
  } else {
    0
  }

  # Feature 5: Range diversity (coefficient of variation of ranges)
  range_diversity <- if (mean(ranges, na.rm = TRUE) > 1e-10) {
    sd(ranges, na.rm = TRUE) / mean(ranges, na.rm = TRUE)
  } else {
    0
  }

  # Feature 6: Skewness diversity (std of skewness values)
  skewness_diversity <- sd(skewness_vals, na.rm = TRUE)

  # Feature 7: Kurtosis diversity (std of kurtosis values)
  kurtosis_diversity <- sd(kurtosis_vals, na.rm = TRUE)

  # Return named list (removed variance_diversity)
  list(
    mean_diversity = mean_diversity,
    acf1_diversity = acf1_diversity,
    complexity_variance = complexity_variance,
    entropy_diversity = entropy_diversity,
    range_diversity = range_diversity,
    skewness_diversity = skewness_diversity,
    kurtosis_diversity = kurtosis_diversity
  )
}

#' Compute Shannon entropy from binned time series
#'
#' Simple histogram-based Shannon entropy estimator
#'
#' @param x Numeric vector
#' @param bins Number of bins for histogram (default 10)
#' @return Shannon entropy estimate
#' @keywords internal
compute_shannon_entropy_simple <- function(x, bins = 10) {
  # Remove NA
  x <- x[!is.na(x)]

  if (length(x) < 10) return(NA)

  # Discretize into bins
  x_breaks <- seq(min(x), max(x), length.out = bins + 1)
  x_binned <- cut(x, breaks = x_breaks, include.lowest = TRUE, labels = FALSE)

  # Compute probabilities
  counts <- table(x_binned)
  probs <- counts / sum(counts)

  # Shannon entropy: H = -sum(p * log(p))
  entropy <- -sum(probs * log(probs + 1e-10))

  return(entropy)
}
