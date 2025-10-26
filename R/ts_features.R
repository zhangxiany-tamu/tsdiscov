#' @useDynLib tsdiscov, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats Box.test HoltWinters acf ar as.ts coef embed lm median na.contiguous na.omit na.pass pacf predict quantile residuals sd stl ts AIC IQR anova ccf fft pchisq pnorm shapiro.test spec.pgram var cor.test ks.test wilcox.test p.adjust fisher.test
#' @importFrom utils head
#' @importFrom graphics hist
NULL

#' Extract all time series features
#'
#' Extracts a set of features from a univariate time series.
#'
#' @param x Numeric vector or ts object representing the time series
#' @param features Character vector of feature sets to extract. Options include
#'   "stats", "acf", "pacf", "trend", "structure", "fft", "scaling", "stl",
#'   "entropy", "canonical", "complexity", and many others. Use "all" (default) for complete set.
#' @return Named list of feature values
#' @export
#' @examples
#' \dontrun{
#' x <- rnorm(100)
#' features <- ts_features(x)
#' }
ts_features <- function(x, features = "all") {
  # Extract frequency information if x is a ts object
  ts_frequency <- NULL
  if (inherits(x, "ts")) {
    ts_frequency <- frequency(x)
    x <- as.numeric(x)
  }

  if (!is.numeric(x)) {
    stop("Input must be a numeric vector or ts object")
  }

  if (length(x) < 3) {
    stop("Time series must have at least 3 observations")
  }

  # Remove NA values
  x <- x[!is.na(x)]

  # Re-check length after NA removal
  if (length(x) < 3) {
    stop("Time series must have at least 3 non-NA observations")
  }

  # Expand "all" to all feature sets from registry
  if (length(features) == 1 && features == "all") {
    features <- get_feature_sets()
  }

  # Validate feature names using registry
  features <- validate_feature_sets(features, strict = TRUE)

  # Get registry
  registry <- .ts_feature_sets_registry()

  # Extract features using registry
  result <- list()
  for (feat_name in features) {
    if (feat_name %in% names(registry)) {
      feat_func <- registry[[feat_name]]
      # Pass frequency to STL features if available
      if (feat_name == "stl" && !is.null(ts_frequency)) {
        result <- c(result, feat_func(x, frequency = ts_frequency))
      } else {
        result <- c(result, feat_func(x))
      }
    }
  }

  result
}

#' Extract Time Series Features as Data Frame
#'
#' Wrapper around ts_features() that returns results as a 1-row data.frame.
#' This format is convenient for extracting features from multiple time series
#' and combining them with rbind().
#'
#' @param x Numeric vector or ts object representing the time series
#' @param features Character vector of feature sets to extract. Use "all" (default)
#'   for the complete set, or specify subsets like c("stats", "acf", "trend").
#' @return A 1-row data.frame with one column per feature, in stable order
#' @export
#' @examples
#' \dontrun{
#' # Extract features from a single series
#' x <- rnorm(100)
#' df <- ts_features_df(x, features = c("stats", "acf"))
#'
#' # Extract from multiple series and combine
#' series_list <- list(rnorm(100), rnorm(100), rnorm(100))
#' results <- do.call(rbind, lapply(series_list, ts_features_df))
#' }
ts_features_df <- function(x, features = "all") {
  result <- ts_features(x, features = features)
  as.data.frame(t(unlist(result)))
}

#' Extract statistical features
#'
#' @param x Numeric vector
#' @return Named list of statistical features
#' @export
ts_stats <- function(x) {
  list(
    mean = cpp_mean(x),
    std = cpp_std(x),
    # variance removed - redundant with std (variance = std^2)
    skewness = cpp_skewness(x),
    kurtosis = cpp_kurtosis(x),
    q25 = cpp_quantile(x, 0.25),
    q50 = cpp_quantile(x, 0.50),
    q75 = cpp_quantile(x, 0.75)
  )
}

#' Extract autocorrelation features
#'
#' @param x Numeric vector
#' @param max_lag Maximum lag for ACF calculation
#' @return Named list of ACF features
#' @export
ts_acf <- function(x, max_lag = 20) {
  acf_vals <- cpp_acf_features(x, max_lag)

  list(
    acf_first_min = acf_vals[1],
    acf_timescale = acf_vals[2],
    acf_sum10 = acf_vals[3],
    acf_lag1 = acf_vals[4],
    acf_lag10 = acf_vals[5],
    x_acf10 = acf_vals[6]
  )
}

#' Extract canonical time series characteristics
#'
#' Computes 20 canonical features covering distribution, autocorrelation,
#' spectral, forecasting, dynamics, and symbolic characteristics.
#'
#' @param x Numeric vector
#' @return Named list of 20 canonical features
#' @export
ts_catch22 <- function(x) {
  list(
    # Distribution features
    mode_5 = cpp_histogram_mode_5(x),
    mode_10 = cpp_histogram_mode_10(x),
    outlier_timing_pos = cpp_outlier_timing_pos(x),
    outlier_timing_neg = cpp_outlier_timing_neg(x),

    # ACF features - removed acf_timescale and acf_first_min (duplicates of ts_acf)

    # Spectral features
    low_freq_power = cpp_welch_power_low_freq(x),
    centroid_freq = cpp_spectral_centroid(x),

    # Forecasting features
    forecast_error = cpp_forecast_error_mean3(x),
    whiten_timescale = cpp_timescale_ratio_after_whitening(x),

    # HRV and dynamics
    high_fluctuation = cpp_high_fluctuation_prop(x),
    trev = cpp_time_reversibility(x),

    # Stretch features
    stretch_high = cpp_stretch_high(x),
    stretch_decreasing = cpp_stretch_low(x),

    # Symbolic/motif features
    entropy_pairs = cpp_motif_three_quantile(x),
    transition_variance = cpp_transition_variance(x, 3),

    # Automutual information
    ami2 = cpp_automutual_info_lag2(x, 5),
    ami_timescale = cpp_automutual_info_first_min(x, 40, 10),

    # Embedding and scaling
    embedding_dist = cpp_embedding_dist_exp_fit(x),
    rs_range = cpp_rs_range(x),
    dfa = cpp_dfa(x),

    # Periodicity
    periodicity = cpp_periodicity_wang(x)
  )
}

#' Extract entropy features
#'
#' @param x Numeric vector
#' @param m Embedding dimension for sample/approximate entropy
#' @param r Tolerance for sample/approximate entropy
#' @return Named list of entropy features
#' @export
ts_entropy <- function(x, m = 2, r = 0.2) {
  list(
    sample_entropy = cpp_sample_entropy(x, m, r),
    approx_entropy = cpp_approximate_entropy(x, m, r),
    perm_entropy = cpp_permutation_entropy(x, 3, 1),
    shannon_entropy = cpp_shannon_entropy(x, 10)
  )
}

#' Extract PACF features
#'
#' @param x Numeric vector
#' @param max_lag Maximum lag for PACF calculation
#' @return Named list of PACF features
#' @export
ts_pacf <- function(x, max_lag = 20) {
  pacf_vals <- cpp_pacf_features(x, max_lag)

  list(
    pacf_lag1 = pacf_vals["pacf_lag1"],
    pacf_lag5 = pacf_vals["pacf_lag5"],
    pacf_lag10 = pacf_vals["pacf_lag10"],
    first_sig_pacf = pacf_vals["first_sig_pacf"],
    sum_sq_pacf = pacf_vals["sum_sq_pacf"],
    x_pacf5 = pacf_vals["x_pacf5"]
  )
}

#' Extract extended PACF features at longer lags
#'
#' Computes partial autocorrelation at longer lags for detecting
#' long-range dependencies and seasonal patterns.
#'
#' @param x Numeric vector
#' @param lags Integer vector of lags to extract (default: c(20, 50))
#' @return Named list with PACF values at specified lags
#' @export
ts_pacf_extended <- function(x, lags = c(20, 50)) {
  n <- length(x)

  # Default result
  default_result <- vector("list", length(lags))
  names(default_result) <- paste0("pacf_lag", lags)
  for (i in 1:length(lags)) default_result[[i]] <- NA_real_

  # Need sufficient data
  if (n < 20) {
    return(default_result)
  }

  # Compute PACF up to max lag
  max_lag <- max(lags)
  if (max_lag >= n / 2) {
    max_lag <- floor(n / 2) - 1
  }

  pacf_result <- tryCatch({
    pacf_full <- cpp_pacf(x, max_lag)

    result <- list()
    for (lag in lags) {
      name <- paste0("pacf_lag", lag)
      if (lag <= length(pacf_full)) {
        result[[name]] <- pacf_full[lag]
      } else {
        result[[name]] <- NA_real_
      }
    }
    result
  }, error = function(e) {
    default_result
  })

  pacf_result
}

#' Extract Permutation Entropy features
#'
#' Computes permutation entropy (Bandt-Pompe) which measures complexity
#' via ordinal patterns. Higher values indicate more randomness.
#'
#' Permutation entropy quantifies the diversity of ordinal patterns in the time series.
#' It's robust to noise and useful for detecting deterministic vs stochastic dynamics.
#'
#' @param x Numeric vector
#' @param dims Integer vector of embedding dimensions (default: c(4, 5))
#'   Note: dim=3 is already provided by ts_entropy() as perm_entropy
#' @param delay Time delay for embedding (default: 1)
#' @return Named list with permutation entropy for each dimension
#' @export
ts_permutation_entropy <- function(x, dims = c(4, 5), delay = 1) {
  n <- length(x)

  # Default result
  default_result <- vector("list", length(dims))
  names(default_result) <- paste0("permutation_entropy_", dims)
  for (i in 1:length(dims)) default_result[[i]] <- NA_real_

  # Need sufficient data
  if (n < max(dims) * delay + 1) {
    return(default_result)
  }

  # Compute permutation entropy for each dimension
  pe_result <- tryCatch({
    result <- list()
    for (dim in dims) {
      name <- paste0("permutation_entropy_", dim)
      result[[name]] <- cpp_permutation_entropy(x, as.integer(dim), as.integer(delay))
    }
    result
  }, error = function(e) {
    default_result
  })

  pe_result
}

#' Extract Multi-scale Entropy features
#'
#' Computes sample entropy at multiple temporal scales via coarse-graining.
#' Reveals complexity across different time scales.
#'
#' At scale k, the series is coarse-grained by averaging non-overlapping windows
#' of size k, then sample entropy is computed on the coarse-grained series.
#'
#' @param x Numeric vector
#' @param scales Integer vector of scales (default: c(2, 5))
#' @param m Embedding dimension for sample entropy (default: 2)
#' @param r Tolerance fraction (default: 0.2)
#' @return Named list with multi-scale entropy at each scale
#' @export
ts_multiscale_entropy <- function(x, scales = c(2, 5), m = 2, r = 0.2) {
  n <- length(x)

  # Default result
  default_result <- vector("list", length(scales))
  names(default_result) <- paste0("multiscale_entropy_", scales)
  for (i in 1:length(scales)) default_result[[i]] <- NA_real_

  # Need sufficient data
  min_n_required <- max(scales) * (m + 10)
  if (n < min_n_required) {
    return(default_result)
  }

  # Compute multi-scale entropy
  mse_result <- tryCatch({
    mse <- cpp_multiscale_entropy(x, as.integer(scales), as.integer(m), r)
    as.list(mse)
  }, error = function(e) {
    default_result
  })

  mse_result
}

#' Extract trend features
#'
#' @param x Numeric vector
#' @return Named list of trend features
#' @export
ts_trend <- function(x) {
  trend_vals <- cpp_linear_trend(x)

  list(
    trend_slope = trend_vals["slope"],
    trend_intercept = trend_vals["intercept"],
    trend_r_squared = trend_vals["r_squared"],
    trend_stderr = trend_vals["stderr"],
    mean_abs_change = cpp_mean_abs_change(x),
    mean_change = cpp_mean_change(x),
    mean_second_deriv = cpp_mean_second_derivative(x)
  )
}

#' Extract time-weighted trend features
#'
#' Computes linear trend with exponential weighting favoring recent data.
#' Useful for detecting recent trend changes while still considering history.
#'
#' @param x Numeric vector
#' @param decay Decay factor for exponential weights (default: 0.95)
#' @return Named list with time-weighted slope and intercept
#' @export
ts_trend_extended <- function(x, decay = 0.95) {
  n <- length(x)

  # Default result
  default_result <- list(
    time_weighted_slope = NA_real_,
    time_weighted_intercept = NA_real_
  )

  if (n < 10) {
    return(default_result)
  }

  # Compute time-weighted trend
  trend_result <- tryCatch({
    cpp_time_weighted_trend(x, decay = decay)
  }, error = function(e) {
    default_result
  })

  as.list(trend_result)
}

#' Extract robust trend features
#'
#' Computes trend using Theil-Sen estimator (median of pairwise slopes).
#' More robust to outliers than ordinary least squares.
#'
#' Theil-Sen is a non-parametric estimator with high breakdown point.
#' It can tolerate up to 29.3% outliers without losing accuracy.
#'
#' @param x Numeric vector
#' @param max_pairs Maximum number of pairs to sample (default: 5000)
#' @return Named list with robust slope, intercept, and slope std
#' @export
ts_robust_trend <- function(x, max_pairs = 5000) {
  n <- length(x)

  # Default result
  default_result <- list(
    robust_slope = NA_real_,
    robust_intercept = NA_real_,
    robust_slope_std = NA_real_
  )

  if (n < 10) {
    return(default_result)
  }

  # Compute robust trend
  trend_result <- tryCatch({
    cpp_robust_trend(x, as.integer(max_pairs))
  }, error = function(e) {
    default_result
  })

  as.list(trend_result)
}

#' Extract Recurrence Quantification Analysis features
#'
#' Analyzes recurrence plots to quantify nonlinear dynamics.
#' Reveals deterministic structure, laminar states, and complexity.
#'
#' RQA is powerful for detecting:
#' - Determinism: presence of predictable dynamics
#' - Laminarity: presence of laminar (slowly changing) states
#' - Complexity: via diagonal line entropy
#'
#' @param x Numeric vector
#' @param threshold_percent Threshold as percentile of distances (default: 0.1)
#' @param min_line_length Minimum line length to count (default: 2)
#' @param max_points Maximum points to analyze (default: 2000)
#' @return Named list with RQA measures
#' @export
ts_recurrence <- function(x, threshold_percent = 0.1, min_line_length = 2, max_points = 2000) {
  n <- length(x)

  # Default result
  default_result <- list(
    recurrence_rate = NA_real_,
    determinism = NA_real_,
    laminarity = NA_real_,
    longest_diagonal = NA_real_,
    entropy_diagonal = NA_real_
  )

  # Need sufficient data
  if (n < 50) {
    return(default_result)
  }

  # Compute RQA
  rqa_result <- tryCatch({
    cpp_recurrence_analysis(x, threshold_percent, as.integer(min_line_length), as.integer(max_points))
  }, error = function(e) {
    default_result
  })

  rqa_result
}

#' Extract structure features
#'
#' @param x Numeric vector
#' @return Named list of structure features
#' @export
ts_structure <- function(x) {
  list(
    # Peak and structure
    number_peaks = cpp_number_peaks(x, 3),
    longest_strike_above = cpp_longest_strike_above_mean(x),
    longest_strike_below = cpp_longest_strike_below_mean(x),

    # Counts
    count_above_mean = cpp_count_above_mean(x),
    count_below_mean = cpp_count_below_mean(x),

    # Ratios and crossings
    ratio_beyond_1sigma = cpp_ratio_beyond_r_sigma(x, 1.0),
    ratio_beyond_2sigma = cpp_ratio_beyond_r_sigma(x, 2.0),
    number_crossings = cpp_number_crossings(x, 0.0),

    # Additional statistics
    absolute_sum_changes = cpp_absolute_sum_of_changes(x),
    range = cpp_range(x),
    mad = cpp_median_absolute_deviation(x),
    coef_variation = cpp_coefficient_of_variation(x),
    benford_correlation = cpp_benford_correlation(x),
    spectral_entropy = cpp_spectral_entropy(x)
  )
}

#' Extract FFT coefficient features
#'
#' @param x Numeric vector
#' @param num_coef Number of FFT coefficients to extract (default 3)
#' @return Named list of FFT features
#' @export
ts_fft <- function(x, num_coef = 3) {
  # Extract individual FFT coefficients (magnitude, real, imaginary, angle)
  fft_coefs <- cpp_fft_coefficients(x, num_coef)

  # Extract aggregated FFT statistics
  fft_agg <- cpp_fft_aggregated(x)

  # Combine into single list
  c(as.list(fft_coefs), as.list(fft_agg))
}

#' Extract scaling and stationarity features
#'
#' @param x Numeric vector
#' @param window_size Window size for stability and lumpiness (default 10)
#' @return Named list of scaling features
#' @export
ts_scaling <- function(x, window_size = 10) {
  list(
    hurst_exponent = cpp_hurst_exponent(x),
    stability = cpp_stability(x, window_size),
    lumpiness = cpp_lumpiness(x, window_size)
  )
}

#' Extract differencing features
#'
#' @param x Numeric vector
#' @param include_ndiffs Logical, whether to include ndiffs (requires urca package)
#' @return Named list of differencing ACF/PACF features
#' @export
ts_diff <- function(x, include_ndiffs = TRUE) {
  diff_features <- cpp_diff_acf_features(x)
  result <- as.list(diff_features)

  # Add ndiffs if requested
  if (include_ndiffs) {
    ndiffs_val <- tryCatch({
      if (!requireNamespace("urca", quietly = TRUE)) {
        NA_real_
      } else {
        # Iteratively difference and test for stationarity using KPSS
        # Returns number of differences needed (0, 1, or 2)
        nd <- NA_real_
        for (d in 0:2) {
          if (d == 0) {
            y <- x
          } else {
            y <- diff(x, differences = d)
          }

          if (length(y) < 10) {
            nd <- NA_real_
            break
          }

          # KPSS test: H0 = stationary
          kpss_result <- urca::ur.kpss(y)
          kpss_stat <- kpss_result@teststat[1]

          # Critical value at 5% significance for 'mu' type is 0.463
          # If stat < critical value, series is stationary
          if (kpss_stat < 0.463) {
            nd <- d
            break
          }
        }
        # If no break occurred (still not stationary after 2 differences), return 2
        if (is.na(nd) && length(y) >= 10) nd <- 2
        nd
      }
    }, error = function(e) NA_real_)

    result <- c(result, list(ndiffs = ndiffs_val))
  }

  return(result)
}

#' Extract STL decomposition features
#'
#' @param x Numeric vector
#' @param frequency Seasonal period (default auto-detected or 1 for non-seasonal)
#' @return Named list of STL features
#' @export
ts_stl <- function(x, frequency = NULL) {
  n <- length(x)

  # Base-R frequency detector (periodogram + ACF cross-check)
  detect_frequency_base <- function(x, min_period = 2, prominence_ratio = 5, acf_min = 0.2) {
    nx <- length(x)
    if (nx < 10) return(1L)

    # Periodogram peak
    sp <- tryCatch(stats::spec.pgram(x, taper = 0.1, detrend = TRUE, plot = FALSE), error = function(e) NULL)
    period_pg <- NA_integer_
    reliable_pg <- FALSE
    if (!is.null(sp) && is.numeric(sp$spec) && is.numeric(sp$freq)) {
      p <- sp$spec; f <- sp$freq
      ok <- is.finite(p) & is.finite(f) & f > 0
      p <- p[ok]; f <- f[ok]
      if (length(p) >= 3) {
        idx <- which.max(p)
        peak <- p[idx]
        medp <- stats::median(p, na.rm = TRUE)
        period_pg <- as.integer(round(1 / f[idx]))
        if (is.finite(period_pg) && period_pg >= min_period && period_pg <= floor(nx / 2)) {
          if (is.finite(peak) && is.finite(medp) && medp > 0 && (peak / medp) >= prominence_ratio) {
            reliable_pg <- TRUE
          }
        }
      }
    }

    # ACF local maximum cross-check
    ac <- tryCatch(stats::acf(x, lag.max = floor(nx / 2), plot = FALSE)$acf, error = function(e) NULL)
    period_acf <- NA_integer_
    reliable_acf <- FALSE
    if (!is.null(ac) && length(ac) > 2) {
      lags <- seq_along(ac) - 1L
      ac <- ac[lags > 0]; lags <- lags[lags > 0]
      locmax <- which(diff(sign(diff(c(-Inf, ac, -Inf)))) == -2)
      if (length(locmax) > 0) {
        cand <- locmax[ac[locmax] >= acf_min & lags[locmax] >= min_period]
        if (length(cand) > 0) {
          cand <- cand[order(-ac[cand], lags[cand])]
          period_acf <- as.integer(lags[cand[1L]])
          reliable_acf <- TRUE
        }
      }
    }

    # Combine
    if (reliable_pg && reliable_acf && is.finite(period_pg) && is.finite(period_acf)) {
      if (abs(period_pg - period_acf) <= 1L) return(as.integer(round((period_pg + period_acf) / 2)))
      return(as.integer(period_pg))
    }
    if (reliable_pg && is.finite(period_pg)) return(as.integer(period_pg))
    if (reliable_acf && is.finite(period_acf)) return(as.integer(period_acf))
    1L
  }

  # Auto-detect frequency if not provided
  if (is.null(frequency)) {
    frequency <- detect_frequency_base(x)
  }

  # If no seasonality or series too short, return NA for all features
  if (frequency == 1 || n < 2 * frequency) {
    return(list(
      trend_strength = NA_real_,
      seasonal_strength = NA_real_,
      spike = NA_real_,
      linearity = NA_real_,
      curvature = NA_real_,
      seas_acf1 = NA_real_,
      seas_pacf = NA_real_,
      e_acf_first_min = NA_real_,
      e_acf_timescale = NA_real_,
      e_acf_sum10 = NA_real_,
      e_acf1 = NA_real_,
      e_acf10 = NA_real_,
      e_diff1_acf1 = NA_real_,
      e_diff1_acf10 = NA_real_,
      e_diff2_acf1 = NA_real_,
      e_diff2_acf10 = NA_real_,
      e_diff1x_pacf5 = NA_real_,
      e_diff2x_pacf5 = NA_real_
    ))
  }

  # Perform STL decomposition
  tryCatch({
    ts_obj <- ts(x, frequency = frequency)
    stl_result <- stl(ts_obj, s.window = "periodic", robust = TRUE)

    trend <- as.numeric(stl_result$time.series[, "trend"])
    seasonal <- as.numeric(stl_result$time.series[, "seasonal"])
    remainder <- as.numeric(stl_result$time.series[, "remainder"])

    # Compute ACF features on remainder
    # cpp_acf_features returns: [acf_first_min, acf_timescale, acf_sum10, acf_lag1, acf_lag10, x_acf10]
    remainder_acf <- cpp_acf_features(remainder, max_lag = 10)
    remainder_diff <- cpp_diff_acf_features(remainder)

    # Compute seasonal ACF/PACF features
    seas_acf1_val <- tryCatch({
      acf_result <- acf(seasonal, lag.max = frequency, plot = FALSE)$acf
      if (length(acf_result) > 1) acf_result[2] else NA_real_  # lag 1
    }, error = function(e) NA_real_)

    seas_pacf_val <- tryCatch({
      pacf_result <- pacf(seasonal, lag.max = frequency, plot = FALSE)$acf
      if (length(pacf_result) >= frequency) pacf_result[frequency] else NA_real_
    }, error = function(e) NA_real_)

    # Combine STL features with remainder ACF features
    list(
      trend_strength = cpp_trend_strength(trend, remainder),
      seasonal_strength = cpp_seasonal_strength(seasonal, remainder),
      spike = cpp_spike(remainder),
      linearity = cpp_linearity(trend),
      curvature = cpp_curvature(trend),
      seas_acf1 = seas_acf1_val,
      seas_pacf = seas_pacf_val,
      # Remainder ACF features (prefix with e_ for STL remainder)
      # Access by index since cpp_acf_features returns unnamed vector
      e_acf_first_min = unname(remainder_acf[1]),
      e_acf_timescale = unname(remainder_acf[2]),
      e_acf_sum10 = unname(remainder_acf[3]),
      e_acf1 = unname(remainder_acf[4]),
      e_acf10 = unname(remainder_acf[6]),  # Sum of squares of first 10 ACF (x_acf10)
      # Remainder differencing features
      e_diff1_acf1 = unname(remainder_diff["diff1_acf1"]),
      e_diff1_acf10 = unname(remainder_diff["diff1_acf10"]),
      e_diff2_acf1 = unname(remainder_diff["diff2_acf1"]),
      e_diff2_acf10 = unname(remainder_diff["diff2_acf10"]),
      e_diff1x_pacf5 = unname(remainder_diff["diff1x_pacf5"]),
      e_diff2x_pacf5 = unname(remainder_diff["diff2x_pacf5"])
    )
  }, error = function(e) {
    # If STL fails, return NA values for all features
    list(
      trend_strength = NA_real_,
      seasonal_strength = NA_real_,
      spike = NA_real_,
      linearity = NA_real_,
      curvature = NA_real_,
      seas_acf1 = NA_real_,
      seas_pacf = NA_real_,
      e_acf_first_min = NA_real_,
      e_acf_timescale = NA_real_,
      e_acf_sum10 = NA_real_,
      e_acf1 = NA_real_,
      e_acf10 = NA_real_,
      e_diff1_acf1 = NA_real_,
      e_diff1_acf10 = NA_real_,
      e_diff2_acf1 = NA_real_,
      e_diff2_acf10 = NA_real_,
      e_diff1x_pacf5 = NA_real_,
      e_diff2x_pacf5 = NA_real_
    )
  })
}

#' Extract heterogeneity features (ARCH/GARCH)
#'
#' Computes measures of volatility clustering and heterogeneity.
#' Pre-whitens using R's ar() function, then computes ARCH effects.
#' Optionally fits GARCH(1,1) using tseries::garch() and computes GARCH effects.
#'
#' @param x Numeric vector
#' @param fit_garch Logical, whether to fit GARCH model (requires tseries package)
#' @return Named list of heterogeneity features
#' @export
ts_heterogeneity <- function(x, fit_garch = TRUE) {
  n <- length(x)

  # Need sufficient data
  if (n < 20) {
    return(list(
      arch_acf = NA_real_,
      garch_acf = NA_real_,
      arch_r2 = NA_real_,
      garch_r2 = NA_real_
    ))
  }

  tryCatch({
    # Pre-whiten using R's built-in ar() function
    ar_fit <- ar(na.contiguous(x), aic = TRUE)
    
    if (is.null(ar_fit$resid) || length(ar_fit$resid) < 15) {
      return(list(
        arch_acf = NA_real_,
        garch_acf = NA_real_,
        arch_r2 = NA_real_,
        garch_r2 = NA_real_
      ))
    }

    # Get pre-whitened series
    whitened <- na.contiguous(ar_fit$resid)
    
    if (length(whitened) < 15) {
      return(list(
        arch_acf = NA_real_,
        garch_acf = NA_real_,
        arch_r2 = NA_real_,
        garch_r2 = NA_real_
      ))
    }

    # Compute ARCH effects using C++ function
    arch_acf_val <- cpp_arch_acf(whitened)
    
    # Compute ARCH R² using R's lm() on squared series
    arch_r2_val <- tryCatch({
      whitened_demean <- whitened - mean(whitened)
      mat <- embed(whitened_demean^2, 13)  # 12 lags + response
      fit <- lm(mat[, 1] ~ mat[, -1])
      r2 <- summary(fit)$r.squared
      if (is.nan(r2)) 1.0 else r2
    }, error = function(e) NA_real_)

    # Initialize GARCH values
    garch_acf_val <- NA_real_
    garch_r2_val <- NA_real_

    # Fit GARCH model if requested
    if (fit_garch && requireNamespace("tseries", quietly = TRUE)) {
      garch_fit <- tryCatch({
        suppressWarnings(tseries::garch(whitened, trace = FALSE))
      }, error = function(e) NULL)

      if (!is.null(garch_fit)) {
        garch_resid <- na.contiguous(residuals(garch_fit))
        
        if (length(garch_resid) >= 15) {
          # Compute GARCH effects
          garch_acf_val <- cpp_arch_acf(garch_resid)
          
          # Compute GARCH R²
          garch_r2_val <- tryCatch({
            mat_g <- embed(garch_resid^2, 13)
            fit_g <- lm(mat_g[, 1] ~ mat_g[, -1])
            r2_g <- summary(fit_g)$r.squared
            if (is.nan(r2_g)) 1.0 else r2_g
          }, error = function(e) NA_real_)
        }
      }
    }

    list(
      arch_acf = arch_acf_val,
      garch_acf = garch_acf_val,
      arch_r2 = arch_r2_val,
      garch_r2 = garch_r2_val
    )
  }, error = function(e) {
    list(
      arch_acf = NA_real_,
      garch_acf = NA_real_,
      arch_r2 = NA_real_,
      garch_r2 = NA_real_
    )
  })
}

#' Extract statistical test features
#'
#' Computes statistical tests for unit roots, nonlinearity, and ARCH effects.
#' Uses existing well-tested R packages: urca for unit root tests and tseries for nonlinearity.
#'
#' @param x Numeric vector
#' @return Named list of statistical test features
#' @export
ts_stattests <- function(x) {
  n <- length(x)

  # Need sufficient data
  if (n < 20) {
    return(list(
      unitroot_kpss = NA_real_,
      unitroot_pp = NA_real_,
      nonlinearity = NA_real_
    ))
  }

  # KPSS unit root test (uses urca package)
  unitroot_kpss_val <- tryCatch({
    if (requireNamespace("urca", quietly = TRUE)) {
      kpss_result <- urca::ur.kpss(x)
      kpss_result@teststat[1]  # Extract test statistic
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)

  # Phillips-Perron unit root test (uses urca package)
  unitroot_pp_val <- tryCatch({
    if (requireNamespace("urca", quietly = TRUE)) {
      pp_result <- urca::ur.pp(x)
      pp_result@teststat[1]  # Extract test statistic
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)

  # Terasvirta nonlinearity test (uses tseries package)
  nonlinearity_val <- tryCatch({
    if (requireNamespace("tseries", quietly = TRUE)) {
      nl_result <- tseries::terasvirta.test(as.ts(x), type = "Chisq")
      # Scaled statistic for comparability across series lengths
      10 * unname(nl_result$stat) / length(x)
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)

  list(
    unitroot_kpss = unitroot_kpss_val,
    unitroot_pp = unitroot_pp_val,
    nonlinearity = nonlinearity_val
  )
}

#' Extract shift detection features
#'
#' Computes features that detect structural breaks and change points in time series.
#' Detects shifts in level (mean), variance, and distribution (KL divergence).
#'
#' @param x Numeric vector
#' @param window_size Window size for rolling statistics (default 10)
#' @return Named list of shift detection features
#' @export
ts_shifts <- function(x, window_size = 10) {
  result <- cpp_shift_detection(x, width = window_size)
  as.list(result)
}

#' Extract miscellaneous features
#'
#' Computes additional useful features including ACF zero crossing,
#' zero proportion, and first derivative statistics.
#'
#' @param x Numeric vector
#' @return Named list of miscellaneous features
#' @export
ts_misc <- function(x) {
  list(
    firstzero_ac = cpp_firstzero_ac(x),
    zero_proportion = cpp_zero_proportion(x),
    std1st_der = cpp_std1st_der(x)
  )
}

#' Computational time series features
#'
#' Calculate computational features for time series characterization.
#' Includes crossing points, embedding analysis, motif detection,
#' and fluctuation analysis.
#'
#' @param x A numeric vector (time series)
#' @return A named list of computational features
#' @export
ts_compengine <- function(x) {
  # Initialize all features to NA for safety
  result_template <- list(
    crossing_points = NA_real_,
    flat_spots = NA_integer_,
    embed2_incircle_1 = NA_real_,
    embed2_incircle_2 = NA_real_,
    motiftwo_entro3 = NA_real_,
    walker_propcross = NA_real_,
    localsimple_mean1 = NA_real_,
    localsimple_lfitac = NA_real_,
    spreadrandomlocal_meantaul_50 = NA_real_,
    spreadrandomlocal_meantaul_ac2 = NA_real_,
    outlierinclude_mdrmd = NA_real_,
    fluctanal_prop_r1 = NA_real_
  )

  if (length(x) < 10) {
    warning("Time series too short (n=", length(x), ") for compengine features. Requires n >= 10. Returning NAs.",
            call. = FALSE)
    return(result_template)
  }

  # Wrap entire computation in tryCatch to ensure we always return 12 features
  tryCatch({
    # C++ implemented features
    crossing_points_val = tryCatch(cpp_crossing_points(x), error = function(e) NA_real_)
    flat_spots_val = tryCatch(cpp_flat_spots(x), error = function(e) NA_integer_)
    embed2_incircle_1_val = tryCatch(cpp_embed2_incircle(x, 1), error = function(e) NA_real_)
    embed2_incircle_2_val = tryCatch(cpp_embed2_incircle(x, 2), error = function(e) NA_real_)
    motiftwo_entro3_val = tryCatch(cpp_motiftwo_entro3(x), error = function(e) NA_real_)
    walker_propcross_val = tryCatch(cpp_walker_propcross(x), error = function(e) NA_real_)

    # Local simple forecasting features
    localsimple_mean1_val = tryCatch({
      # Use mean of past 1 value to predict next
      n = length(x)
      res = numeric(n - 1)
      for (i in 2:n) {
        res[i-1] = x[i-1] - x[i]  # Prediction (x[i-1]) - actual (x[i])
      }
      # First zero crossing of ACF of residuals
      cpp_firstzero_ac(res)
    }, error = function(e) NA_real_)

    localsimple_lfitac_val = tryCatch({
      # Use local linear fit to predict next value
      # tau = first zero crossing of ACF
      tau = cpp_firstzero_ac(x)
      if (is.na(tau) || tau == 0) tau = 3  # Default to 3 if no zero crossing

      n = length(x)
      if (tau >= n - 1) {
        NA_real_  # Don't use return() here - it exits the entire function!
      } else {
        evalr = (tau + 1):n
        res = numeric(length(evalr))

        for (i in seq_along(evalr)) {
          idx = evalr[i]
          window = x[(idx - tau):(idx - 1)]
          time_idx = 1:tau

          # Check if window has sufficient variation for linear fit
          if (sd(window, na.rm = TRUE) < .Machine$double.eps) {
            # Constant window - use mean as prediction
            res[i] = mean(window, na.rm = TRUE) - x[idx]
          } else {
            # Fit linear model
            fit = lm(window ~ time_idx)
            # Check if model is well-defined
            if (length(coef(fit)) == 2 && !any(is.na(coef(fit)))) {
              pred = predict(fit, newdata = data.frame(time_idx = tau + 1))
              res[i] = pred - x[idx]
            } else {
              # Fall back to mean prediction
              res[i] = mean(window, na.rm = TRUE) - x[idx]
            }
          }
        }

        # First zero crossing of ACF of residuals
        cpp_firstzero_ac(res)
      }  # End of else block
    }, error = function(e) NA_real_)

    # Bootstrap stationarity features
    spreadrandomlocal_meantaul_50_val = tryCatch({
      spreadrandomlocal_meantaul(x, 50)
    }, error = function(e) NA_real_, warning = function(w) {
      NA_real_
    })

    spreadrandomlocal_meantaul_ac2_val = tryCatch({
      tau_ac2 = 2 * cpp_firstzero_ac(x)
      if (is.na(tau_ac2) || tau_ac2 == 0) {
        NA_real_  # Don't use return() here - it exits the entire function!
      } else {
        spreadrandomlocal_meantaul(x, tau_ac2)
      }
    }, error = function(e) NA_real_, warning = function(w) NA_real_)

    # Outlier inclusion median
    outlierinclude_mdrmd_val = tryCatch({
      outlierinclude_mdrmd(x)
    }, error = function(e) NA_real_)

    # Fluctuation analysis
    fluctanal_prop_r1_val = tryCatch({
      fluctanal_prop_r1(x)
    }, error = function(e) NA_real_)

    # Build result list - ensure all 12 features are present
    result <- list(
      crossing_points = crossing_points_val,
      flat_spots = flat_spots_val,
      embed2_incircle_1 = embed2_incircle_1_val,
      embed2_incircle_2 = embed2_incircle_2_val,
      motiftwo_entro3 = motiftwo_entro3_val,
      walker_propcross = walker_propcross_val,
      localsimple_mean1 = localsimple_mean1_val,
      localsimple_lfitac = localsimple_lfitac_val,
      spreadrandomlocal_meantaul_50 = spreadrandomlocal_meantaul_50_val,
      spreadrandomlocal_meantaul_ac2 = spreadrandomlocal_meantaul_ac2_val,
      outlierinclude_mdrmd = outlierinclude_mdrmd_val,
      fluctanal_prop_r1 = fluctanal_prop_r1_val
    )

    # Validate result has correct structure
    if (!is.list(result) || length(result) != 12) {
      warning("ts_compengine: Unexpected result structure (length=", length(result),
              "). Returning template with NAs.", call. = FALSE)
      return(result_template)
    }

    return(result)

  }, error = function(e) {
    warning("ts_compengine: Error during feature extraction: ", e$message,
            ". Returning NAs for all features.", call. = FALSE)
    return(result_template)
  })
}

# Helper function: Bootstrap stationarity measure
# Note: This function uses random sampling and will produce slightly
# different results on each run (~2% variation). Set seed for reproducibility.
spreadrandomlocal_meantaul <- function(y, l) {
  numSegs = 100
  n = length(y)

  # Check if series is long enough
  min_required_length = ceiling(l / 0.9)
  if (l > 0.9 * n) {
    warning("Time series too short (n=", n, ") for spreadrandomlocal_meantaul with window l=", l,
            ". Requires n >= ", min_required_length, ". Returning NA.",
            call. = FALSE)
    return(NA_real_)
  }

  qs = numeric(numSegs)
  for (j in 1:numSegs) {
    # Pick random start point
    ist = sample.int(n - l, 1)
    ifh = ist + l - 1
    ysub = y[ist:ifh]

    # First zero crossing of ACF
    taul = cpp_firstzero_ac(ysub)
    qs[j] = taul
  }

  mean(qs, na.rm = TRUE)
}

# Helper function: Outlier inclusion median
outlierinclude_mdrmd <- function(y, zscored = TRUE) {
  if (length(unique(y)) == 1) {
    return(NA_real_)
  }

  # Z-score if requested
  if (zscored) {
    y = scale(y)[,1]
    isd = 1
  } else {
    isd = sd(y, na.rm = TRUE)
  }

  n = length(y)
  inc = 0.01 * isd
  thr = seq(from = 0, to = max(abs(y), na.rm = TRUE), by = inc)

  if (length(thr) == 0) return(NA_real_)

  msDt = numeric(length(thr))
  msDtp = numeric(length(thr))

  for (i in seq_along(thr)) {
    th = thr[i]
    r = which(abs(y) >= th)

    if (length(r) < 2) {
      msDt[i] = NA
      msDtp[i] = 0
      next
    }

    Dt_exc = diff(r)
    msDt[i] = median(r) / (n / 2) - 1
    msDtp[i] = length(Dt_exc) / n * 100
  }

  # Trim where less than 2% of data included
  trimthr = 2
  valid_idx = which(msDtp > trimthr)
  if (length(valid_idx) == 0) return(NA_real_)

  mj = max(valid_idx)
  msDt = msDt[1:mj]

  median(msDt, na.rm = TRUE)
}

# Helper function: Fluctuation analysis
fluctanal_prop_r1 <- function(x) {
  q = 2
  tauStep = 50
  k = 1
  n = length(x)

  # Replace NA with 0
  x_NA0 = ifelse(!is.na(x), x, 0)
  y = cumsum(x_NA0)

  # Exponentially spaced scales
  taur = unique(round(exp(seq(from = log(5), to = log(floor(n / 2)), length.out = tauStep))))
  ntau = length(taur)

  if (ntau < 8) {
    return(NA_real_)
  }

  Fl = numeric(ntau)

  for (i in 1:ntau) {
    tau = taur[i]

    # Buffer time series at scale tau
    y_buff = split(y, ceiling(seq_along(y) / tau))

    # Remove trailing partial segment
    if (length(y_buff[[length(y_buff)]]) < tau) {
      y_buff = y_buff[-length(y_buff)]
    }

    if (length(y_buff) == 0) {
      Fl[i] = NA
      next
    }

    tt = 1:tau

    # Detrend each segment
    y_dt = numeric(length(y_buff))
    for (j in seq_along(y_buff)) {
      if (length(y_buff[[j]]) == tau) {
        fit = lm(y_buff[[j]] ~ tt)
        resid = residuals(fit)
        y_dt[j] = max(resid) - min(resid)
      } else {
        y_dt[j] = NA
      }
    }

    # Fluctuation function
    Fl[i] = (mean(y_dt^q, na.rm = TRUE))^(1 / q)
  }

  logtt = log(taur)
  logFF = log(Fl)
  ntt = ntau

  # Find breakpoint between two scaling regimes
  sserr = rep(NA, ntt)
  minPoints = 6

  for (i in minPoints:(ntt - minPoints)) {
    r1 = 1:i
    p1 = lm(logFF[r1] ~ logtt[r1])

    r2 = i:ntt
    p2 = lm(logFF[r2] ~ logtt[r2])

    sserr[i] = sum(residuals(p1)^2) + sum(residuals(p2)^2)
  }

  breakPt = which.min(sserr)
  r1 = 1:breakPt

  # Proportion in first regime
  prop_r1 = length(r1) / ntt
  return(prop_r1)
}

#' Exponential smoothing parameters
#'
#' Estimate smoothing parameters using base R's HoltWinters function
#'
#' @param x A numeric vector (time series)
#' @param frequency Seasonal period (default auto-detected or 1 for non-seasonal)
#' @return A named list of Holt and Holt-Winters parameters
#' @export
ts_holt <- function(x, frequency = NULL) {
  n <- length(x)

  # Auto-detect frequency if not provided
  if (is.null(frequency)) {
    if (n >= 24) {
      frequency <- 12  # Monthly default
    } else if (n >= 14) {
      frequency <- 7   # Weekly default
    } else {
      frequency <- 1   # No seasonality
    }
  }

  # Holt parameters (no seasonal component)
  holt_alpha <- NA_real_
  holt_beta <- NA_real_

  if (n >= 10) {
    holt_result <- tryCatch({
      # Check data suitability
      if (any(is.na(x)) || sd(x, na.rm = TRUE) < .Machine$double.eps) {
        return(list(alpha = NA_real_, beta = NA_real_))
      }

      ts_obj <- ts(x, frequency = frequency)
      # Suppress expected warnings for white noise/unsuitable data
      fit <- suppressWarnings(
        HoltWinters(ts_obj, gamma = FALSE, optim.control = list(maxit = 100))
      )

      if (is.null(fit) || any(is.na(c(fit$alpha, fit$beta)))) {
        list(alpha = NA_real_, beta = NA_real_)
      } else {
        list(alpha = fit$alpha, beta = fit$beta)
      }
    }, error = function(e) list(alpha = NA_real_, beta = NA_real_))

    holt_alpha <- holt_result$alpha
    holt_beta <- holt_result$beta
  }

  # Holt-Winters parameters (with seasonal component)
  hw_alpha <- NA_real_
  hw_beta <- NA_real_
  hw_gamma <- NA_real_

  if (frequency > 1 && n >= 2 * frequency) {
    hw_result <- tryCatch({
      ts_obj <- ts(x, frequency = frequency)

      # Check data suitability for HoltWinters
      if (any(is.na(x)) || sd(x, na.rm = TRUE) < .Machine$double.eps) {
        return(list(alpha = NA_real_, beta = NA_real_, gamma = NA_real_))
      }

      # Try additive seasonal model first
      # Suppress optimization warnings (not errors) as they indicate
      # the data doesn't have clear seasonal patterns
      fit <- tryCatch({
        suppressWarnings(
          HoltWinters(ts_obj, seasonal = "additive", optim.control = list(maxit = 100))
        )
      }, error = function(e) {
        # If additive fails, try multiplicative (if all values are positive)
        if (all(x > 0)) {
          tryCatch({
            suppressWarnings(
              HoltWinters(ts_obj, seasonal = "multiplicative", optim.control = list(maxit = 100))
            )
          }, error = function(e2) NULL)
        } else {
          NULL
        }
      })

      if (is.null(fit)) {
        list(alpha = NA_real_, beta = NA_real_, gamma = NA_real_)
      } else {
        list(alpha = fit$alpha, beta = fit$beta, gamma = fit$gamma)
      }
    }, error = function(e) list(alpha = NA_real_, beta = NA_real_, gamma = NA_real_))

    hw_alpha <- hw_result$alpha
    hw_beta <- hw_result$beta
    hw_gamma <- hw_result$gamma
  }

  list(
    holt_alpha = holt_alpha,
    holt_beta = holt_beta,
    hw_alpha = hw_alpha,
    hw_beta = hw_beta,
    hw_gamma = hw_gamma
  )
}

#' Extract complexity and nonlinearity features
#'
#' Computes advanced complexity measures including C3 nonlinearity,
#' CID complexity, Lempel-Ziv complexity, and temporal distribution features.
#'
#' @param x Numeric vector
#' @return Named list of complexity features
#' @export
ts_tsfresh <- function(x) {
  list(
    # C3 non-linearity (lag 1)
    c3_lag1 = cpp_c3(x, lag = 1),

    # CID complexity
    cid_ce = cpp_cid_ce(x, normalize = TRUE),

    # Lempel-Ziv complexity
    lempel_ziv = cpp_lempel_ziv_complexity(x, bins = 10),

    # Index mass quantiles
    index_mass_q25 = cpp_index_mass_quantile(x, q = 0.25),
    index_mass_q50 = cpp_index_mass_quantile(x, q = 0.50),
    index_mass_q75 = cpp_index_mass_quantile(x, q = 0.75),

    # Change quantiles (middle 50%)
    change_quantiles = cpp_change_quantiles(x, ql = 0.25, qh = 0.75, isabs = TRUE, f_agg = "mean"),

    # Mean of top 3 absolute maxima
    mean_abs_max_3 = cpp_mean_n_absolute_max(x, number_of_maxima = 3),

    # Energy ratio in first chunk (out of 10)
    energy_ratio_first = cpp_energy_ratio_by_chunks(x, num_segments = 10, segment_focus = 0),

    # Fourier entropy
    fourier_entropy = cpp_fourier_entropy(x, bins = 10)
  )
}

#' Extract supplementary features
#'
#' Computes additional features including aggregations,
#' uniqueness measures, and distribution characteristics.
#'
#' @param x Numeric vector
#' @return Named list of supplementary features
#' @export
ts_tsfresh_supp <- function(x) {
  list(
    # Simple aggregations
    abs_energy = cpp_abs_energy(x),
    sum_values = cpp_sum_values(x),

    # Uniqueness and duplicates
    has_duplicate = cpp_has_duplicate(x),
    has_duplicate_max = cpp_has_duplicate_max(x),
    has_duplicate_min = cpp_has_duplicate_min(x),
    percentage_reoccurring = cpp_percentage_reoccurring(x),
    sum_reoccurring = cpp_sum_reoccurring(x),
    ratio_unique_values = cpp_ratio_unique_values(x),

    # Distribution characteristics (boolean checks)
    symmetry_looking = cpp_symmetry_looking(x, r = 0.05),
    large_standard_deviation = cpp_large_standard_deviation(x, r = 0.25),
    variance_larger_than_std = cpp_variance_larger_than_std(x)
  )
}

#' Extract Continuous Wavelet Transform features
#'
#' Computes continuous wavelet transform using Ricker (Mexican hat) wavelet
#' for multi-scale time-frequency analysis. Uses direct implementation without
#' requiring external wavelet packages.
#'
#' @param x Numeric vector
#' @param widths Widths for CWT (default c(2, 5, 10, 20))
#' @return Named list of CWT features
#' @export
ts_cwt <- function(x, widths = c(2, 5, 10, 20)) {
  n <- length(x)

  if (n < 10) {
    # Return NA values for too-short series
    result <- list()
    for (w in widths) {
      result[[paste0("cwt_coeff_width", w)]] <- NA_real_
    }
    return(result)
  }

  # Perform CWT using Ricker (Mexican hat) wavelet
  # Implementation uses direct convolution without external packages
  result <- list()

  tryCatch({
    for (w in widths) {
      # Create Ricker wavelet at this width
      # Use continuous wavelet transform
      # For simplicity, we'll compute the maximum CWT coefficient at this scale

      # Approximate Ricker wavelet CWT
      # Ricker wavelet: (1 - x^2/a^2) * exp(-x^2/(2a^2))
      half_width <- as.integer(3 * w)
      t_vals <- seq(-half_width, half_width)

      # Ricker wavelet
      ricker <- (2 / (sqrt(3 * w) * pi^0.25)) * (1 - (t_vals^2 / w^2)) * exp(-t_vals^2 / (2 * w^2))

      # Convolve with series
      if (length(x) >= length(ricker)) {
        conv_result <- stats::convolve(x, rev(ricker), type = "open")
        # Take max absolute coefficient
        max_coeff <- max(abs(conv_result), na.rm = TRUE)
        result[[paste0("cwt_coeff_width", w)]] <- max_coeff
      } else {
        result[[paste0("cwt_coeff_width", w)]] <- NA_real_
      }
    }
  }, error = function(e) {
    for (w in widths) {
      result[[paste0("cwt_coeff_width", w)]] <- NA_real_
    }
  })

  return(result)
}

#' Extract basic statistical features
#'
#' Computes basic statistics: maximum, minimum, and root mean square.
#'
#' @param x Numeric vector
#' @return Named list of basic statistics
#' @export
ts_basic <- function(x) {
  list(
    maximum = max(x, na.rm = TRUE),
    minimum = min(x, na.rm = TRUE),
    root_mean_square = sqrt(mean(x^2, na.rm = TRUE))
  )
}

#' Extract extrema location features
#'
#' Computes normalized locations (0 to 1) of first and last occurrences
#' of maximum and minimum values.
#'
#' @param x Numeric vector
#' @return Named list of extrema locations
#' @export
ts_extrema <- function(x) {
  n <- length(x)

  if (n == 0) {
    return(list(
      first_loc_max = NA_real_,
      last_loc_max = NA_real_,
      first_loc_min = NA_real_,
      last_loc_min = NA_real_
    ))
  }

  # Remove NAs for finding extrema
  x_clean <- x[!is.na(x)]
  n_clean <- length(x_clean)

  if (n_clean == 0) {
    return(list(
      first_loc_max = NA_real_,
      last_loc_max = NA_real_,
      first_loc_min = NA_real_,
      last_loc_min = NA_real_
    ))
  }

  # Find locations (normalized to 0-1)
  list(
    first_loc_max = which.max(x_clean) / n_clean,
    last_loc_max = (n_clean - which.max(rev(x_clean)) + 1) / n_clean,
    first_loc_min = which.min(x_clean) / n_clean,
    last_loc_min = (n_clean - which.min(rev(x_clean)) + 1) / n_clean
  )
}

#' Extract AR (AutoRegressive) coefficient features
#'
#' Fits AR models at different orders and extracts coefficients.
#' Uses base R stats::ar() function with Yule-Walker estimation.
#'
#' @param x Numeric vector
#' @param orders Vector of AR orders to extract (default: c(1, 2, 5, 10))
#' @return Named list of AR coefficients
#' @export
ts_ar <- function(x, orders = c(1, 2, 5, 10)) {
  n <- length(x)

  # Initialize result with NAs
  result <- list()
  for (p in orders) {
    result[[paste0("ar_coef_", p)]] <- NA_real_
  }

  # Need enough data for AR fitting
  if (n < 10) {
    return(result)
  }

  # Check for sufficient variation
  if (sd(x, na.rm = TRUE) < .Machine$double.eps) {
    return(result)
  }

  # Fit AR model
  ar_result <- tryCatch({
    # Use Yule-Walker method (default in ar())
    # Set order.max to highest requested order
    ar_fit <- ar(x, aic = FALSE, order.max = max(orders), na.action = na.pass)

    # Extract coefficients at requested orders
    for (p in orders) {
      if (p <= length(ar_fit$ar)) {
        result[[paste0("ar_coef_", p)]] <- ar_fit$ar[p]
      }
      # else leave as NA (already initialized)
    }
    result
  }, error = function(e) {
    result  # Return NAs on error
  })

  ar_result
}

#' Extract ADF (Augmented Dickey-Fuller) test features
#'
#' Performs Augmented Dickey-Fuller unit root test using internal implementation.
#' Tests null hypothesis: time series has a unit root (non-stationary).
#'
#' The test uses automatic lag selection based on Schwert (1989):
#' lag = floor(12 * (n/100)^(1/4))
#'
#' @param x Numeric vector
#' @param max_lag Maximum lag order. If -1 (default), automatically selected.
#' @param include_trend Include trend in test regression (default: TRUE)
#' @return Named list with ADF test statistic and p-value
#' @export
ts_adf <- function(x, max_lag = -1, include_trend = TRUE) {
  n <- length(x)

  # Default return for edge cases
  default_result <- list(adf_stat = NA_real_, adf_pvalue = NA_real_)

  # Need sufficient data
  if (n < 10) {
    return(default_result)
  }

  # Remove NAs
  x <- x[!is.na(x)]
  if (length(x) < 10) {
    return(default_result)
  }

  # Perform ADF test using Rcpp implementation
  adf_result <- tryCatch({
    test <- cpp_adf_test(x, max_lag = max_lag, include_trend = include_trend)
    list(
      adf_stat = test$statistic,
      adf_pvalue = test$p.value
    )
  }, error = function(e) {
    default_result
  })

  adf_result
}

#' Extract Welch Power Spectral Density features
#'
#' Computes Welch's power spectral density estimate using overlapping windows.
#' More robust than simple FFT for noisy signals. Returns power in frequency bins.
#'
#' Uses Hanning window with 50% overlap by default for optimal variance reduction.
#'
#' @param x Numeric vector
#' @param n_freq Number of frequency bins (default: 5)
#' @param window_size Window size for segments. If -1, automatically determined.
#' @param overlap Overlap fraction between windows (default: 0.5)
#' @return Named list with power in each frequency bin
#' @export
ts_welch <- function(x, n_freq = 5, window_size = -1, overlap = 0.5) {
  n <- length(x)

  # Default return for edge cases
  default_result <- vector("list", n_freq)
  names(default_result) <- paste0("welch_psd_bin", 1:n_freq)
  for (i in 1:n_freq) default_result[[i]] <- NA_real_

  # Need sufficient data
  if (n < 10) {
    return(default_result)
  }

  # Remove NAs
  x <- x[!is.na(x)]
  if (length(x) < 10) {
    return(default_result)
  }

  # Compute Welch PSD using Rcpp implementation
  welch_result <- tryCatch({
    psd <- cpp_welch_psd(x, n_freq = n_freq, window_size = window_size, overlap = overlap)
    as.list(psd)
  }, error = function(e) {
    default_result
  })

  welch_result
}

#' Extract Matrix Profile features
#'
#' Computes simplified matrix profile for pattern discovery.
#' Matrix profile finds nearest neighbor distances for all subsequences.
#' Low matrix profile values indicate motifs (repeated patterns).
#' High matrix profile values indicate discords (anomalies).
#'
#' @param x Numeric vector
#' @param window_size Subsequence length. If -1, automatically determined (10\% of series length).
#' @return Named list with matrix profile statistics and pattern scores
#' @export
ts_matrix_profile <- function(x, window_size = -1) {
  n <- length(x)

  # Default return for edge cases
  default_result <- list(
    matrix_profile_min = NA_real_,
    matrix_profile_mean = NA_real_,
    matrix_profile_max = NA_real_,
    motif_score = NA_real_,
    discord_score = NA_real_
  )

  # Need sufficient data
  if (n < 20) {
    return(default_result)
  }

  # Remove NAs
  if (any(is.na(x))) {
    return(default_result)
  }

  # Compute Matrix Profile using Rcpp implementation
  mp_result <- tryCatch({
    cpp_matrix_profile(x, window_size = window_size)
  }, error = function(e) {
    default_result
  })

  mp_result
}

#' Extract Friedrich coefficient features
#'
#' Computes Friedrich coefficients from Langevin model fitting.
#' Models drift as polynomial: dx/dt = a₀ + a₁·x + a₂·x² + a₃·x³.
#'
#' These coefficients capture deterministic dynamics in noisy time series.
#' Used in physics, climate science, and neuroscience for detecting
#' nonlinear drift terms.
#'
#' @param x Numeric vector
#' @param max_order Maximum polynomial order (default: 3)
#' @return Named list with polynomial coefficients a₀, a₁, a₂, a₃
#' @export
ts_friedrich <- function(x, max_order = 3) {
  n <- length(x)

  # Default return for edge cases
  default_result <- vector("list", max_order + 1)
  names(default_result) <- paste0("friedrich_coef_", 0:max_order)
  for (i in 1:(max_order + 1)) default_result[[i]] <- NA_real_

  # Need sufficient data
  if (n < 20) {
    return(default_result)
  }

  # Remove NAs
  if (any(is.na(x))) {
    return(default_result)
  }

  # Compute Friedrich coefficients using Rcpp implementation
  friedrich_result <- tryCatch({
    cpp_friedrich_coefficients(x, max_order = max_order)
  }, error = function(e) {
    default_result
  })

  friedrich_result
}

#' Extract Langevin fixed point features
#'
#' Finds fixed points (stable states) from the Langevin model.
#' Fixed points are where the drift function (from Friedrich coefficients) equals zero.
#'
#' These identify equilibrium states in dynamical systems.
#' Used in physics and neuroscience for stability analysis.
#'
#' @param x Numeric vector
#' @param max_order Maximum polynomial order for Friedrich model (default: 3)
#' @return Named list with first fixed point and maximum fixed point
#' @export
ts_langevin <- function(x, max_order = 3) {
  n <- length(x)

  # Default return
  default_result <- list(
    langevin_fixed_point = NA_real_,
    langevin_max_fixed_point = NA_real_
  )

  # Need sufficient data
  if (n < 20) {
    return(default_result)
  }

  # Remove NAs
  if (any(is.na(x))) {
    return(default_result)
  }

  # Compute fixed points using Rcpp implementation
  langevin_result <- tryCatch({
    cpp_langevin_fixed_point(x, max_order = max_order)
  }, error = function(e) {
    default_result
  })

  langevin_result
}

#' Ljung-Box test for serial correlation
#'
#' Tests for autocorrelation in residuals using the Ljung-Box test statistic.
#' Returns test statistics and p-values at lags 10 and 20.
#'
#' @param x Numeric vector
#' @return Named list with Ljung-Box statistics and p-values
#' @export
ts_ljungbox <- function(x) {
  n <- length(x)

  # Default return
  default_result <- list(
    ljungbox_stat_10 = NA_real_,
    ljungbox_pval_10 = NA_real_,
    ljungbox_stat_20 = NA_real_,
    ljungbox_pval_20 = NA_real_
  )

  # Need sufficient data
  if (n < 25) {
    return(default_result)
  }

  # Test at lag 10
  test10 <- tryCatch({
    Box.test(x, lag = 10, type = "Ljung-Box")
  }, error = function(e) NULL)

  # Test at lag 20
  test20 <- tryCatch({
    Box.test(x, lag = 20, type = "Ljung-Box")
  }, error = function(e) NULL)

  list(
    ljungbox_stat_10 = if (!is.null(test10)) unname(test10$statistic) else NA_real_,
    ljungbox_pval_10 = if (!is.null(test10)) test10$p.value else NA_real_,
    ljungbox_stat_20 = if (!is.null(test20)) unname(test20$statistic) else NA_real_,
    ljungbox_pval_20 = if (!is.null(test20)) test20$p.value else NA_real_
  )
}

#' Spectral shape descriptors
#'
#' Computes shape characteristics of the power spectral density including
#' flatness, slope, rolloff, and bandwidth.
#'
#' @param x Numeric vector
#' @return Named list with spectral shape features
#' @export
ts_spectral_shape <- function(x) {
  n <- length(x)

  # Default return
  default_result <- list(
    spectral_flatness = NA_real_,
    spectral_slope = NA_real_,
    spectral_rolloff_95 = NA_real_,
    spectral_bandwidth = NA_real_
  )

  # Need sufficient data
  if (n < 20) {
    return(default_result)
  }

  # Remove NAs
  x <- x[!is.na(x)]
  if (length(x) < 20) {
    return(default_result)
  }

  # Compute spectral shape using Rcpp
  shape_result <- tryCatch({
    cpp_spectral_shape(x)
  }, error = function(e) {
    default_result
  })

  shape_result
}

#' Poincaré plot features
#'
#' Computes SD1, SD2, and their ratio from lag-1 embedding.
#'
#' @param x Numeric vector
#' @return Named list with poincare_sd1, poincare_sd2, poincare_ratio
#' @export
ts_poincare <- function(x) {
  n <- length(x)

  default_result <- list(
    poincare_sd1 = NA_real_,
    poincare_sd2 = NA_real_,
    poincare_ratio = NA_real_
  )

  if (n < 10) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 10) return(default_result)

  tryCatch({
    x1 <- x[-length(x)]
    x2 <- x[-1]

    sd1 <- sqrt(var(x2 - x1) / 2)
    sd2 <- sqrt(var(x2 + x1) / 2)
    ratio <- if (sd2 > 1e-10) sd1 / sd2 else NA_real_

    list(
      poincare_sd1 = sd1,
      poincare_sd2 = sd2,
      poincare_ratio = ratio
    )
  }, error = function(e) default_result)
}

#' Run-length and up-down statistics
#'
#' Computes statistics based on runs above/below median.
#'
#' @param x Numeric vector
#' @return Named list with run statistics
#' @export
ts_runs <- function(x) {
  n <- length(x)

  default_result <- list(
    runs_n = NA_real_,
    runs_mean_length = NA_real_,
    runs_updown_ratio = NA_real_,
    runs_switch_rate = NA_real_
  )

  if (n < 10) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 10) return(default_result)

  tryCatch({
    m <- median(x, na.rm = TRUE)
    labs <- ifelse(x >= m, 1L, 0L)
    rl <- rle(labs)

    num_runs <- length(rl$lengths)
    mean_run_len <- mean(rl$lengths)

    n_up <- sum(rl$values == 1)
    n_down <- sum(rl$values == 0)
    up_down_ratio <- if (n_down > 0) n_up / n_down else NA_real_

    sw <- sum(diff(labs) != 0, na.rm = TRUE) / (length(labs) - 1)

    list(
      runs_n = num_runs,
      runs_mean_length = mean_run_len,
      runs_updown_ratio = up_down_ratio,
      runs_switch_rate = sw
    )
  }, error = function(e) default_result)
}

#' Wave and impulse shape factors
#'
#' Computes crest, impulse, shape, and clearance factors.
#'
#' @param x Numeric vector
#' @return Named list with shape factors
#' @export
ts_shape_factors <- function(x) {
  n <- length(x)

  default_result <- list(
    shape_crest = NA_real_,
    shape_impulse = NA_real_,
    shape_shape = NA_real_,
    shape_clearance = NA_real_
  )

  if (n < 10) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 10) return(default_result)

  tryCatch({
    xa <- abs(x)
    RMS <- sqrt(mean(xa^2))
    A <- mean(xa)
    R <- mean(sqrt(xa))
    P <- max(xa)

    crest <- if (RMS > 1e-10) P / RMS else NA_real_
    impulse <- if (A > 1e-10) P / A else NA_real_
    shape <- if (A > 1e-10) RMS / A else NA_real_
    clearance <- if (R > 1e-10) P / (R^2) else NA_real_

    list(
      shape_crest = crest,
      shape_impulse = impulse,
      shape_shape = shape,
      shape_clearance = clearance
    )
  }, error = function(e) default_result)
}

#' Forecastability measure
#'
#' Computes 1 - normalized spectral entropy.
#'
#' @param x Numeric vector
#' @return Named list with forecastability
#' @export
ts_forecastability <- function(x) {
  n <- length(x)

  default_result <- list(forecastability = NA_real_)

  if (n < 20) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 20) return(default_result)

  tryCatch({
    sp <- spec.pgram(x, taper = 0, detrend = TRUE, plot = FALSE)
    p <- sp$spec
    p <- p / sum(p)

    H <- -sum(p[p > 0] * log(p[p > 0]))
    Hn <- H / log(length(p))
    forecastability <- 1 - Hn

    list(forecastability = forecastability)
  }, error = function(e) default_result)
}

#' Robust frequency detection
#'
#' Detects dominant frequency using combined periodogram and ACF analysis.
#' Uses prominence ratio and ACF threshold to ensure reliable detection.
#'
#' @param x Numeric vector
#' @param min_period Minimum period to consider (default: 2)
#' @param max_period Maximum period to consider (default: floor(n/2))
#' @param prominence_ratio Peak must exceed median by this factor (default: 5)
#' @param acf_min Minimum ACF value for peak detection (default: 0.2)
#' @return Named list with freq_est (estimated period) and nperiods (number of cycles)
#' @export
ts_frequency <- function(x, min_period = 2, max_period = NULL,
                        prominence_ratio = 5, acf_min = 0.2) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)

  if (n < 10) return(list(freq_est = 1L, nperiods = 0L))

  if (is.null(max_period)) max_period <- floor(n / 2)
  max_period <- max(min_period, min(max_period, floor(n / 2)))

  # 1) Periodogram-based period (primary)
  sp <- tryCatch(spec.pgram(x, taper = 0.1, fast = TRUE, detrend = TRUE,
                            plot = FALSE), error = function(e) NULL)
  period_est_pg <- NA_integer_
  reliable_pg <- FALSE

  if (!is.null(sp) && is.numeric(sp$spec) && is.numeric(sp$freq)) {
    p <- sp$spec
    f <- sp$freq  # cycles per sample, in (0, 0.5]

    # Exclude DC and invalids
    valid <- is.finite(p) & is.finite(f) & f > 0
    p <- p[valid]
    f <- f[valid]

    if (length(p) >= 3) {
      # Pick dominant peak
      idx <- which.max(p)
      peak <- p[idx]
      medp <- median(p, na.rm = TRUE)

      # Convert to period in samples
      period_pg <- round(1 / f[idx])

      # Clamp to allowed range
      if (is.finite(period_pg) && period_pg >= min_period && period_pg <= max_period) {
        # Prominence check: dominant peak should stand out vs median
        if (is.finite(peak) && is.finite(medp) && medp > 0 && (peak / medp) >= prominence_ratio) {
          period_est_pg <- as.integer(period_pg)
          reliable_pg <- TRUE
        }
      }
    }
  }

  # 2) ACF-based period (secondary / cross-check)
  period_est_acf <- NA_integer_
  reliable_acf <- FALSE

  ac <- tryCatch(acf(x, lag.max = floor(n / 2), plot = FALSE)$acf,
                 error = function(e) NULL)

  if (!is.null(ac) && length(ac) > 2) {
    lags <- seq_along(ac) - 1L
    ac <- ac[lags > 0]
    lags <- lags[lags > 0]

    # Find local maxima
    locmax <- which(diff(sign(diff(c(-Inf, ac, -Inf)))) == -2)

    if (length(locmax) > 0) {
      # Choose strongest candidate above acf_min and within allowable periods
      cand <- locmax[ac[locmax] >= acf_min & lags[locmax] >= min_period &
                     lags[locmax] <= max_period]

      if (length(cand) > 0) {
        # Among candidates, take highest ACF, then smallest lag
        cand <- cand[order(-ac[cand], lags[cand])]
        period_est_acf <- as.integer(lags[cand[1L]])
        reliable_acf <- TRUE
      }
    }
  }

  # 3) Combine estimates
  if (reliable_pg && reliable_acf && is.finite(period_est_pg) && is.finite(period_est_acf)) {
    # If both reliable and close → average
    if (abs(period_est_pg - period_est_acf) <= 1L) {
      period <- as.integer(round((period_est_pg + period_est_acf) / 2))
    } else {
      # Prefer periodogram if much more prominent
      period <- as.integer(period_est_pg)
    }
  } else if (reliable_pg) {
    period <- as.integer(period_est_pg)
  } else if (reliable_acf) {
    period <- as.integer(period_est_acf)
  } else {
    period <- 1L
  }

  if (!is.finite(period) || period < min_period || period > max_period) {
    period <- 1L
  }

  nperiods <- if (period > 1L) as.integer(floor(n / period)) else 0L

  list(freq_est = period, nperiods = nperiods)
}

#' ARIMA residual diagnostics
#'
#' Fits simple ARIMA, tests residuals for whiteness, normality, and ARCH effects.
#'
#' @param x Numeric vector
#' @return Named list with diagnostic p-values
#' @export
ts_arima_diag <- function(x) {
  n <- length(x)

  default_result <- list(
    arima_ljungbox_pval = NA_real_,
    arima_normality_pval = NA_real_,
    arima_arch_pval = NA_real_
  )

  if (n < 30) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 30) return(default_result)

  tryCatch({
    # Fit simple ARIMA(1,0,1)
    fit <- suppressWarnings(stats::arima(x, order = c(1, 0, 1)))
    res <- residuals(fit)

    # Ljung-Box test
    lag <- min(10, floor(length(res) / 5))
    lb <- Box.test(res, type = "Ljung-Box", lag = lag, fitdf = 2)
    lb_pval <- lb$p.value

    # Normality test
    norm_pval <- if (length(res) <= 5000) {
      shapiro.test(res)$p.value
    } else {
      s <- mean(scale(res)^3)
      k <- mean(scale(res)^4)
      JB <- length(res) / 6 * (s^2 + (k - 3)^2 / 4)
      1 - pchisq(JB, df = 2)
    }

    # ARCH LM test
    y <- res^2
    K <- min(5, floor(length(res) / 10))
    if (length(y) > K + 1) {
      mat <- embed(y, K + 1)
      r2 <- summary(lm(mat[, 1] ~ mat[, -1]))$r.squared
      LM <- nrow(mat) * r2
      arch_pval <- 1 - pchisq(LM, K)
    } else {
      arch_pval <- NA_real_
    }

    list(
      arima_ljungbox_pval = lb_pval,
      arima_normality_pval = norm_pval,
      arima_arch_pval = arch_pval
    )
  }, error = function(e) default_result)
}

#' Hjorth parameters
#'
#' Computes Hjorth activity, mobility, and complexity parameters.
#'
#' @param x Numeric vector
#' @return Named list with activity, mobility, complexity
#' @export
ts_hjorth <- function(x) {
  n <- length(x)

  default_result <- list(
    hjorth_activity = NA_real_,
    hjorth_mobility = NA_real_,
    hjorth_complexity = NA_real_
  )

  if (n < 10) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 10) return(default_result)

  tryCatch({
    activity <- var(x)
    dx <- diff(x)
    var_dx <- var(dx)
    mobility <- sqrt(var_dx / activity)

    ddx <- diff(dx)
    var_ddx <- var(ddx)
    complexity <- sqrt(var_ddx / var_dx) / mobility

    list(
      hjorth_activity = activity,
      hjorth_mobility = mobility,
      hjorth_complexity = complexity
    )
  }, error = function(e) default_result)
}

#' Line length
#'
#' Sum of absolute differences (total variation).
#'
#' @param x Numeric vector
#' @return Named list with line_length
#' @export
ts_line_length <- function(x) {
  n <- length(x)

  default_result <- list(line_length = NA_real_)

  if (n < 2) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 2) return(default_result)

  tryCatch({
    list(line_length = sum(abs(diff(x))))
  }, error = function(e) default_result)
}

#' Slope sign changes
#'
#' Number of slope sign changes in the series.
#'
#' @param x Numeric vector
#' @return Named list with slope_sign_changes
#' @export
ts_ssc <- function(x) {
  n <- length(x)

  default_result <- list(slope_sign_changes = NA_real_)

  if (n < 3) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 3) return(default_result)

  tryCatch({
    dx <- diff(x)
    # Remove zeros before sign computation
    dx_nonzero <- dx[dx != 0]
    if (length(dx_nonzero) < 2) {
      return(list(slope_sign_changes = 0))
    }

    sign_changes <- sum(diff(sign(dx_nonzero)) != 0)

    list(slope_sign_changes = sign_changes)
  }, error = function(e) default_result)
}

#' Turning points
#'
#' Rate of local extrema in the series.
#'
#' @param x Numeric vector
#' @return Named list with turning_point_rate
#' @export
ts_turning_points <- function(x) {
  n <- length(x)

  default_result <- list(turning_point_rate = NA_real_)

  if (n < 3) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 3) return(default_result)

  tryCatch({
    count <- 0
    for (i in 2:(length(x) - 1)) {
      if ((x[i] - x[i-1]) * (x[i+1] - x[i]) < 0) {
        count <- count + 1
      }
    }

    rate <- count / (length(x) - 2)

    list(turning_point_rate = rate)
  }, error = function(e) default_result)
}

#' Monotonicity
#'
#' Kendall correlation with time index and Mann-Kendall p-value.
#'
#' @param x Numeric vector
#' @return Named list with kendall_tau_abs, mk_pvalue
#' @export
ts_monotonicity <- function(x) {
  n <- length(x)

  default_result <- list(
    kendall_tau_abs = NA_real_,
    mk_pvalue = NA_real_
  )

  if (n < 10) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 10) return(default_result)

  tryCatch({
    test <- cor.test(x, seq_along(x), method = "kendall")

    list(
      kendall_tau_abs = abs(test$estimate),
      mk_pvalue = test$p.value
    )
  }, error = function(e) default_result)
}

#' Quartile coefficient of dispersion
#'
#' Normalized measure of spread based on quartiles.
#'
#' @param x Numeric vector
#' @return Named list with quartile_coef_dispersion
#' @export
ts_qcd <- function(x) {
  n <- length(x)

  default_result <- list(quartile_coef_dispersion = NA_real_)

  if (n < 4) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 4) return(default_result)

  tryCatch({
    q <- quantile(x, c(0.25, 0.75), type = 7)
    qcd <- (q[2] - q[1]) / (q[2] + q[1])

    list(quartile_coef_dispersion = qcd)
  }, error = function(e) default_result)
}

#' IQR-based properties
#'
#' IQR relative to median and outlier proportion.
#'
#' @param x Numeric vector
#' @return Named list with iqr_over_median_abs, outlier_prop_iqr
#' @export
ts_iqr_props <- function(x) {
  n <- length(x)

  default_result <- list(
    iqr_over_median_abs = NA_real_,
    outlier_prop_iqr = NA_real_
  )

  if (n < 4) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 4) return(default_result)

  tryCatch({
    q <- quantile(x, c(0.25, 0.5, 0.75), type = 7)
    iqr <- q[3] - q[1]

    iqr_ratio <- iqr / (abs(q[2]) + 1e-10)

    lower <- q[1] - 1.5 * iqr
    upper <- q[3] + 1.5 * iqr
    outlier_prop <- mean(x < lower | x > upper)

    list(
      iqr_over_median_abs = iqr_ratio,
      outlier_prop_iqr = outlier_prop
    )
  }, error = function(e) default_result)
}

#' Gini coefficient
#'
#' Measure of statistical dispersion.
#'
#' @param x Numeric vector
#' @return Named list with gini_coefficient
#' @export
ts_gini <- function(x) {
  n <- length(x)

  default_result <- list(gini_coefficient = NA_real_)

  if (n < 2) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 2) return(default_result)

  tryCatch({
    # Use absolute values for stability
    x_abs <- abs(x)
    n <- length(x_abs)

    # Compute Gini: sum of all pairwise absolute differences
    total_diff <- 0
    for (i in 1:n) {
      for (j in 1:n) {
        total_diff <- total_diff + abs(x_abs[i] - x_abs[j])
      }
    }

    gini <- total_diff / (2 * n * sum(x_abs))

    list(gini_coefficient = gini)
  }, error = function(e) default_result)
}

#' Average magnitude difference function
#'
#' AMDF-based features for periodicity detection.
#'
#' @param x Numeric vector
#' @param max_lag Maximum lag to consider
#' @return Named list with amdf_min_lag, amdf_min_value, amdf_ratio
#' @export
ts_amdf <- function(x, max_lag = NULL) {
  n <- length(x)

  default_result <- list(
    amdf_min_lag = NA_real_,
    amdf_min_value = NA_real_,
    amdf_ratio = NA_real_
  )

  if (n < 20) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 20) return(default_result)

  tryCatch({
    if (is.null(max_lag)) max_lag <- floor(length(x) / 2)
    max_lag <- min(max_lag, floor(length(x) / 2))

    amdf <- vapply(1:max_lag, function(k) {
      mean(abs(x[(k+1):length(x)] - x[1:(length(x)-k)]))
    }, 0.0)

    min_idx <- which.min(amdf)
    min_val <- amdf[min_idx]
    med_amdf <- median(amdf)

    list(
      amdf_min_lag = min_idx,
      amdf_min_value = min_val,
      amdf_ratio = if (med_amdf > 0) min_val / med_amdf else NA_real_
    )
  }, error = function(e) default_result)
}

#' ACF integral and peak features
#'
#' ACF-based integral timescale and peak prominence.
#'
#' @param x Numeric vector
#' @param max_lag Maximum lag
#' @return Named list with acf_integral_tau, acf_peak_prominence
#' @export
ts_acf_integrals <- function(x, max_lag = NULL) {
  n <- length(x)

  default_result <- list(
    acf_integral_tau = NA_real_,
    acf_peak_prominence = NA_real_
  )

  if (n < 20) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 20) return(default_result)

  tryCatch({
    if (is.null(max_lag)) max_lag <- floor(length(x) / 2)

    ac <- acf(x, lag.max = max_lag, plot = FALSE)$acf[-1]

    # Integral tau: sum positive ACF until first zero crossing
    integral_tau <- 0
    for (i in seq_along(ac)) {
      if (ac[i] <= 0) break
      integral_tau <- integral_tau + ac[i]
    }

    # Peak prominence
    ac_abs <- abs(ac)
    peak <- max(ac[ac > 0], na.rm = TRUE)
    med_ac <- median(ac_abs, na.rm = TRUE)
    prominence <- if (med_ac > 0) peak / med_ac else NA_real_

    list(
      acf_integral_tau = integral_tau,
      acf_peak_prominence = prominence
    )
  }, error = function(e) default_result)
}

#' Spectral moments
#'
#' Higher-order moments of power spectral density.
#'
#' @param x Numeric vector
#' @param n_freq Number of frequencies
#' @return Named list with spectral_spread, spectral_skewness, spectral_kurtosis, dc_ratio
#' @export
ts_spectral_moments <- function(x, n_freq = 128) {
  n <- length(x)

  default_result <- list(
    spectral_spread = NA_real_,
    spectral_skewness = NA_real_,
    spectral_kurtosis = NA_real_,
    dc_ratio = NA_real_
  )

  if (n < 20) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 20) return(default_result)

  tryCatch({
    sp <- spec.pgram(x, taper = 0.1, detrend = TRUE, plot = FALSE)
    p <- sp$spec / sum(sp$spec)  # Normalize
    f <- sp$freq

    # Spectral centroid
    mu <- sum(f * p)

    # Higher moments
    spread <- sqrt(sum((f - mu)^2 * p))
    skewness <- if (spread > 1e-10) sum(((f - mu)^3) * p) / spread^3 else NA_real_
    kurtosis <- if (spread > 1e-10) sum(((f - mu)^4) * p) / spread^4 else NA_real_

    # DC ratio (lowest frequency bin / total)
    dc_ratio <- p[1]

    list(
      spectral_spread = spread,
      spectral_skewness = skewness,
      spectral_kurtosis = kurtosis,
      dc_ratio = dc_ratio
    )
  }, error = function(e) default_result)
}

#' Peak-to-sidelobe ratio
#'
#' Ratio of dominant peak to second largest peak in spectrum.
#'
#' @param x Numeric vector
#' @return Named list with peak_sidelobe_ratio
#' @export
ts_spectral_psr <- function(x) {
  n <- length(x)

  default_result <- list(peak_sidelobe_ratio = NA_real_)

  if (n < 20) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 20) return(default_result)

  tryCatch({
    sp <- spec.pgram(x, taper = 0.1, detrend = TRUE, plot = FALSE)
    p <- sp$spec

    # Exclude DC (first bin)
    p_no_dc <- p[-1]

    if (length(p_no_dc) < 2) return(default_result)

    # Find two largest peaks
    sorted <- sort(p_no_dc, decreasing = TRUE)
    psr <- sorted[1] / sorted[2]

    list(peak_sidelobe_ratio = psr)
  }, error = function(e) default_result)
}

#' Teager-Kaiser energy operator
#'
#' Mean and std of Teager-Kaiser energy.
#'
#' @param x Numeric vector
#' @return Named list with teager_mean, teager_std
#' @export
ts_tkao <- function(x) {
  n <- length(x)

  default_result <- list(
    teager_mean = NA_real_,
    teager_std = NA_real_
  )

  if (n < 3) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 3) return(default_result)

  tryCatch({
    # Teager-Kaiser: psi[t] = x[t]^2 - x[t-1]*x[t+1]
    psi <- x[2:(length(x)-1)]^2 - x[1:(length(x)-2)] * x[3:length(x)]

    list(
      teager_mean = mean(psi),
      teager_std = sd(psi)
    )
  }, error = function(e) default_result)
}

#' Time reversal asymmetry at multiple lags
#'
#' Measures asymmetry under time reversal.
#'
#' @param x Numeric vector
#' @param lags Lags to compute (default: c(1, 2, 3))
#' @return Named list with trav_asym_l1, trav_asym_l2, trav_asym_l3
#' @export
ts_time_reversal_multi <- function(x, lags = c(1, 2, 3)) {
  n <- length(x)

  default_result <- list(
    trav_asym_l1 = NA_real_,
    trav_asym_l2 = NA_real_,
    trav_asym_l3 = NA_real_
  )

  if (n < 10) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 10) return(default_result)

  tryCatch({
    sd_x <- sd(x)
    if (sd_x < 1e-10) return(default_result)

    result <- list()
    for (i in seq_along(lags)) {
      tau <- lags[i]
      if (length(x) > tau) {
        x_t <- x[1:(length(x) - tau)]
        x_t_tau <- x[(tau + 1):length(x)]

        asym <- mean((x_t_tau - x_t)^3) / sd_x^3
        result[[paste0("trav_asym_l", tau)]] <- asym
      } else {
        result[[paste0("trav_asym_l", tau)]] <- NA_real_
      }
    }

    result
  }, error = function(e) default_result)
}

#' Simple STL seasonality features
#'
#' Phase and amplitude of seasonal component.
#'
#' @param x Numeric vector
#' @param frequency Seasonal period (default: auto-detect)
#' @return Named list with seasonal_phase, seasonal_amplitude
#' @export
ts_stl_simple <- function(x, frequency = NULL) {
  n <- length(x)

  default_result <- list(
    seasonal_phase = NA_real_,
    seasonal_amplitude = NA_real_
  )

  if (n < 20) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 20) return(default_result)

  tryCatch({
    # Auto-detect frequency if not provided
    if (is.null(frequency)) {
      freq_result <- ts_frequency(x)
      frequency <- freq_result$freq_est
      if (is.na(frequency) || frequency <= 1) frequency <- 1
    }

    # Need at least 2 periods
    if (frequency <= 1 || length(x) < 2 * frequency) {
      return(default_result)
    }

    # Run STL
    ts_obj <- ts(x, frequency = frequency)
    stl_result <- suppressWarnings(stl(ts_obj, s.window = "periodic", robust = TRUE))
    seasonal <- stl_result$time.series[, "seasonal"]

    # Phase: location of maximum
    phase <- which.max(seasonal) / length(seasonal)

    # Amplitude: range
    amplitude <- max(seasonal) - min(seasonal)

    list(
      seasonal_phase = phase,
      seasonal_amplitude = amplitude
    )
  }, error = function(e) default_result)
}

#' CUSUM features
#'
#' Cumulative sum statistics for change detection.
#'
#' @param x Numeric vector
#' @return Named list with cusum_max, cusum_min, cusum_range
#' @export
ts_cusum <- function(x) {
  n <- length(x)

  default_result <- list(
    cusum_max = NA_real_,
    cusum_min = NA_real_,
    cusum_range = NA_real_
  )

  if (n < 2) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 2) return(default_result)

  tryCatch({
    c <- cumsum(x - mean(x))

    list(
      cusum_max = max(c),
      cusum_min = min(c),
      cusum_range = max(c) - min(c)
    )
  }, error = function(e) default_result)
}

#' Spectral Q-factor
#'
#' Quality factor of dominant spectral peak.
#'
#' @param x Numeric vector
#' @return Named list with spectral_q
#' @export
ts_spectral_qfactor <- function(x) {
  n <- length(x)

  default_result <- list(spectral_q = NA_real_)

  if (n < 20) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 20) return(default_result)

  tryCatch({
    sp <- spec.pgram(x, taper = 0.1, detrend = TRUE, plot = FALSE)
    P <- sp$spec
    F <- sp$freq

    i <- which.max(P)
    p0 <- P[i]
    half <- p0 / 2

    # Find crossings around peak
    l <- if (i > 1) max(1, which.max(P[1:i] < half)) else 1
    r <- if (i < length(P)) (i - 1) + which.max(P[i:length(P)] < half) else length(P)

    fw <- if (l > 0 && r > i) F[r] - F[l] else NA_real_
    q_factor <- if (is.finite(fw) && fw > 0) F[i] / fw else NA_real_

    list(spectral_q = q_factor)
  }, error = function(e) default_result)
}

#' Spectral band ratios
#'
#' Power in low, mid, and high frequency bands.
#'
#' @param x Numeric vector
#' @return Named list with band_low_ratio, band_mid_ratio, band_high_ratio
#' @export
ts_spectral_bands <- function(x) {
  n <- length(x)

  default_result <- list(
    band_low_ratio = NA_real_,
    band_mid_ratio = NA_real_,
    band_high_ratio = NA_real_
  )

  if (n < 20) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 20) return(default_result)

  tryCatch({
    sp <- spec.pgram(x, taper = 0.1, detrend = TRUE, plot = FALSE)
    P <- sp$spec
    F <- sp$freq
    Pn <- P / sum(P)

    low <- sum(Pn[F <= 0.1])
    mid <- sum(Pn[F > 0.1 & F <= 0.3])
    high <- sum(Pn[F > 0.3])

    list(
      band_low_ratio = low,
      band_mid_ratio = mid,
      band_high_ratio = high
    )
  }, error = function(e) default_result)
}

#' Spectral crest factor
#'
#' Ratio of peak to mean spectral power.
#'
#' @param x Numeric vector
#' @return Named list with spectral_crest_factor
#' @export
ts_spectral_crest <- function(x) {
  n <- length(x)

  default_result <- list(spectral_crest_factor = NA_real_)

  if (n < 20) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 20) return(default_result)

  tryCatch({
    sp <- spec.pgram(x, taper = 0.1, detrend = TRUE, plot = FALSE)
    P <- sp$spec

    list(spectral_crest_factor = max(P, na.rm = TRUE) / mean(P, na.rm = TRUE))
  }, error = function(e) default_result)
}

#' Hill tail index
#'
#' Estimates tail heaviness via Hill estimator.
#'
#' @param x Numeric vector
#' @return Named list with hill_tail_index
#' @export
ts_tail_index <- function(x) {
  n <- length(x)

  default_result <- list(hill_tail_index = NA_real_)

  if (n < 20) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 20) return(default_result)

  tryCatch({
    y <- abs(x)
    y <- y[y > 0]

    if (length(y) < 20) return(default_result)

    thr <- quantile(y, 0.95, type = 7)
    z <- sort(y[y >= thr], decreasing = TRUE)
    k <- length(z) - 1

    est <- if (k >= 5) 1 / mean(log(z[1:k]) - log(z[k + 1])) else NA_real_

    list(hill_tail_index = est)
  }, error = function(e) default_result)
}

#' Robust shape statistics
#'
#' Bowley skewness and Moors kurtosis based on quantiles.
#'
#' @param x Numeric vector
#' @return Named list with bowley_skew, moors_kurtosis
#' @export
ts_robust_shape <- function(x) {
  n <- length(x)

  default_result <- list(
    bowley_skew = NA_real_,
    moors_kurtosis = NA_real_
  )

  if (n < 8) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 8) return(default_result)

  tryCatch({
    q <- quantile(x, c(0.25, 0.5, 0.75), type = 7)
    bow <- if (q[3] > q[1]) (q[3] + q[1] - 2 * q[2]) / (q[3] - q[1]) else NA_real_

    oct <- quantile(x, seq(0.125, 0.875, by = 0.125), type = 7)
    den <- (oct[6] - oct[2])
    moor <- if (den > 0) (oct[7] - oct[5] + oct[3] - oct[1]) / den else NA_real_

    list(
      bowley_skew = bow,
      moors_kurtosis = moor
    )
  }, error = function(e) default_result)
}

#' ACF summary statistics
#'
#' Mean absolute ACF and decay index.
#'
#' @param x Numeric vector
#' @return Named list with acf_mean_abs_1to10, acf_decay_index
#' @export
ts_acf_summary <- function(x) {
  n <- length(x)

  default_result <- list(
    acf_mean_abs_1to10 = NA_real_,
    acf_decay_index = NA_real_
  )

  if (n < 20) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 20) return(default_result)

  tryCatch({
    ac <- acf(x, lag.max = 10, plot = FALSE)$acf[-1]
    mean_abs <- mean(abs(ac), na.rm = TRUE)

    decay <- which(abs(ac) < 0.1)[1]
    decay <- if (is.na(decay)) NA_integer_ else as.integer(decay)

    list(
      acf_mean_abs_1to10 = mean_abs,
      acf_decay_index = decay
    )
  }, error = function(e) default_result)
}

#' Volatility clustering
#'
#' ACF of absolute and squared values.
#'
#' @param x Numeric vector
#' @return Named list with ac_abs_l1, ac_sq_l1, ac_sq_l5, ac_sq_l10
#' @export
ts_volatility_cluster <- function(x) {
  n <- length(x)

  default_result <- list(
    ac_abs_l1 = NA_real_,
    ac_sq_l1 = NA_real_,
    ac_sq_l5 = NA_real_,
    ac_sq_l10 = NA_real_
  )

  if (n < 20) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 20) return(default_result)

  tryCatch({
    ac_abs <- acf(abs(x), lag.max = 10, plot = FALSE)$acf
    ac_sq <- acf(x^2, lag.max = 10, plot = FALSE)$acf

    get_lag <- function(ac, k) {
      if (length(ac) > k) ac[k + 1] else NA_real_
    }

    list(
      ac_abs_l1 = get_lag(ac_abs, 1),
      ac_sq_l1 = get_lag(ac_sq, 1),
      ac_sq_l5 = get_lag(ac_sq, 5),
      ac_sq_l10 = get_lag(ac_sq, 10)
    )
  }, error = function(e) default_result)
}

#' Peak interval statistics
#'
#' Count and spacing of peaks.
#'
#' @param x Numeric vector
#' @return Named list with peak_count, peak_mean_interval, peak_cv_interval
#' @export
ts_peak_intervals <- function(x) {
  n <- length(x)

  default_result <- list(
    peak_count = 0L,
    peak_mean_interval = NA_real_,
    peak_cv_interval = NA_real_
  )

  if (n < 5) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 5) return(default_result)

  tryCatch({
    # Simple peak detection: strict greater than neighbors
    prom <- median(abs(x - median(x)))
    pk_idx <- c()

    for (i in 2:(length(x) - 1)) {
      if (x[i] > x[i - 1] && x[i] > x[i + 1]) {
        # Check prominence
        if ((x[i] - min(x[i - 1], x[i + 1])) >= 0.1 * prom) {
          pk_idx <- c(pk_idx, i)
        }
      }
    }

    cnt <- length(pk_idx)

    if (cnt >= 2) {
      ints <- diff(pk_idx)
      list(
        peak_count = as.integer(cnt),
        peak_mean_interval = mean(ints),
        peak_cv_interval = sd(ints) / mean(ints)
      )
    } else {
      list(
        peak_count = as.integer(cnt),
        peak_mean_interval = NA_real_,
        peak_cv_interval = NA_real_
      )
    }
  }, error = function(e) default_result)
}

#' Dwell times
#'
#' Mean run lengths above and below median.
#'
#' @param x Numeric vector
#' @return Named list with dwell_above_mean, dwell_below_mean
#' @export
ts_dwell_times <- function(x) {
  n <- length(x)

  default_result <- list(
    dwell_above_mean = NA_real_,
    dwell_below_mean = NA_real_
  )

  if (n < 2) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 2) return(default_result)

  tryCatch({
    lab <- ifelse(x >= median(x, na.rm = TRUE), 1L, 0L)
    rl <- rle(lab)

    above <- rl$lengths[rl$values == 1]
    below <- rl$lengths[rl$values == 0]

    list(
      dwell_above_mean = if (length(above) > 0) mean(above) else NA_real_,
      dwell_below_mean = if (length(below) > 0) mean(below) else NA_real_
    )
  }, error = function(e) default_result)
}

#' Tail ratio statistics
#'
#' Ratios of extreme quantiles.
#'
#' @param x Numeric vector
#' @return Named list with q99_q95, q95_q90
#' @export
ts_tail_ratios <- function(x) {
  n <- length(x)

  default_result <- list(
    q99_q95 = NA_real_,
    q95_q90 = NA_real_
  )

  if (n < 20) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 20) return(default_result)

  tryCatch({
    qs <- quantile(x, c(0.90, 0.95, 0.99), type = 7)

    list(
      q99_q95 = (qs[3] - qs[2]) / (abs(qs[2]) + 1e-12),
      q95_q90 = (qs[2] - qs[1]) / (abs(qs[1]) + 1e-12)
    )
  }, error = function(e) default_result)
}

#' Winsorized statistics
#'
#' Mean and variance after Winsorizing at 5%.
#'
#' @param x Numeric vector
#' @return Named list with winsor_mean, winsor_var
#' @export
ts_winsor_stats <- function(x) {
  n <- length(x)

  default_result <- list(
    winsor_mean = NA_real_,
    winsor_var = NA_real_
  )

  if (n < 10) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 10) return(default_result)

  tryCatch({
    p <- 0.05
    qs <- quantile(x, c(p, 1 - p), type = 7)
    xw <- pmin(pmax(x, qs[1]), qs[2])

    list(
      winsor_mean = mean(xw),
      winsor_var = var(xw)
    )
  }, error = function(e) default_result)
}

#' Rényi entropy (order 2)
#'
#' Rényi entropy of order 2 from histogram.
#'
#' @param x Numeric vector
#' @return Named list with renyi2
#' @export
ts_renyi2_entropy <- function(x) {
  n <- length(x)

  default_result <- list(renyi2 = NA_real_)

  if (n < 10) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 10) return(default_result)

  tryCatch({
    br <- max(8, min(64, floor(sqrt(n))))
    h <- hist(x, breaks = br, plot = FALSE)
    p <- h$density * diff(h$breaks)
    p <- p[p > 0]
    p <- p / sum(p)

    list(renyi2 = -log(sum(p^2)))
  }, error = function(e) default_result)
}

#' Fourier harmonic features
#'
#' Amplitude and phase of first harmonic.
#'
#' @param x Numeric vector
#' @return Named list with harmonic_amp1, harmonic_phase1
#' @export
ts_fourier1 <- function(x) {
  n <- length(x)

  default_result <- list(
    harmonic_amp1 = NA_real_,
    harmonic_phase1 = NA_real_
  )

  if (n < 20) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 20) return(default_result)

  tryCatch({
    # Get period estimate
    freq_result <- ts_frequency(x)
    P <- max(2L, as.integer(freq_result$freq_est))

    if (P < 2 || n < 2 * P) return(default_result)

    t <- seq_len(length(x))
    fit <- lm(x ~ I(sin(2 * pi * t / P)) + I(cos(2 * pi * t / P)))
    coeffs <- coef(fit)

    a <- coeffs[2]
    b <- coeffs[3]

    list(
      harmonic_amp1 = sqrt(a^2 + b^2),
      harmonic_phase1 = atan2(b, a)
    )
  }, error = function(e) default_result)
}

#' AR(1) half-life
#'
#' AR(1) coefficient and implied half-life.
#'
#' @param x Numeric vector
#' @return Named list with ar1_phi, half_life
#' @export
ts_ar1_halflife <- function(x) {
  n <- length(x)

  default_result <- list(
    ar1_phi = NA_real_,
    half_life = NA_real_
  )

  if (n < 10) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 10) return(default_result)

  tryCatch({
    phi <- coef(lm(x[-1] ~ x[-n] - 1))[1]
    phi <- unname(phi)

    hl <- if (is.finite(phi) && abs(phi) < 1) -log(2) / log(abs(phi)) else NA_real_

    list(
      ar1_phi = phi,
      half_life = hl
    )
  }, error = function(e) default_result)
}

#' CUSUM of squared deviations
#'
#' Range of cumulative sum of squared deviations.
#'
#' @param x Numeric vector
#' @return Named list with cusum_sq_range
#' @export
ts_cusum_sq <- function(x) {
  n <- length(x)

  default_result <- list(cusum_sq_range = NA_real_)

  if (n < 2) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 2) return(default_result)

  tryCatch({
    xs <- x - mean(x)
    c2 <- cumsum(xs^2)

    list(cusum_sq_range = max(c2) - min(c2))
  }, error = function(e) default_result)
}

#' Katz fractal dimension
#'
#' Measures complexity via polyline distance metric.
#'
#' @param x Numeric vector
#' @return Named list with katz_fd
#' @export
ts_fractal_katz <- function(x) {
  n <- length(x)

  default_result <- list(katz_fd = NA_real_)

  if (n < 3) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 3) return(default_result)

  tryCatch({
    # Total line length
    L <- sum(abs(diff(x)))

    # Max distance from first point along polyline
    dists <- sqrt(cumsum(diff(x)^2) + (seq_along(diff(x)))^2)
    d <- max(dists)

    katz <- if (d > 0 && L > 0) log(n) / (log(n) + log(d / L)) else NA_real_

    list(katz_fd = katz)
  }, error = function(e) default_result)
}

#' Petrosian fractal dimension
#'
#' Fractal dimension based on sign changes in differences.
#'
#' @param x Numeric vector
#' @return Named list with petrosian_fd
#' @export
ts_fractal_petrosian <- function(x) {
  n <- length(x)

  default_result <- list(petrosian_fd = NA_real_)

  if (n < 3) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 3) return(default_result)

  tryCatch({
    # Number of sign changes in diff(x)
    dx <- diff(x)
    n_delta <- sum(dx[-1] * dx[-length(dx)] < 0)

    petro <- log(n) / (log(n) + log(n / (n + 0.4 * n_delta)))

    list(petrosian_fd = petro)
  }, error = function(e) default_result)
}

#' Higuchi fractal dimension
#'
#' Fractal dimension via average curve length scaling.
#'
#' @param x Numeric vector
#' @return Named list with higuchi_fd
#' @export
ts_fractal_higuchi <- function(x) {
  n <- length(x)

  default_result <- list(higuchi_fd = NA_real_)

  if (n < 20) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 20) return(default_result)

  tryCatch({
    K <- min(10, floor(n / 4))
    if (K < 2) return(default_result)

    L_k <- numeric(K)

    for (k in 1:K) {
      L_m <- numeric(k)
      for (m in 1:k) {
        idx <- seq(m, n, by = k)
        if (length(idx) < 2) next
        L_m[m] <- sum(abs(diff(x[idx]))) * (n - 1) / ((length(idx) - 1) * k)
      }
      L_k[k] <- mean(L_m[L_m > 0])
    }

    valid <- L_k > 0
    if (sum(valid) < 3) return(default_result)

    fit <- lm(log(L_k[valid]) ~ log(1 / (1:K)[valid]))
    higuchi <- -coef(fit)[2]

    list(higuchi_fd = unname(higuchi))
  }, error = function(e) default_result)
}

#' First-order entropy rate
#'
#' Conditional entropy H(X|X-1) from transition matrix.
#'
#' @param x Numeric vector
#' @return Named list with entropy_rate1
#' @export
ts_entropy_rate_order1 <- function(x) {
  n <- length(x)

  default_result <- list(entropy_rate1 = NA_real_)

  if (n < 10) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 10) return(default_result)

  tryCatch({
    B <- 10
    bins <- cut(x, breaks = B, labels = FALSE, include.lowest = TRUE)
    if (any(is.na(bins))) return(default_result)

    # Transition matrix
    trans <- table(bins[-n], bins[-1])
    trans <- trans / sum(trans)

    # H(X|X-1) = sum_i,j p(i,j) log(p(j|i))
    marginal <- rowSums(trans)
    conditional <- trans / (marginal + 1e-10)

    rate <- 0
    for (i in 1:nrow(trans)) {
      for (j in 1:ncol(trans)) {
        if (trans[i, j] > 0) {
          rate <- rate - trans[i, j] * log(conditional[i, j])
        }
      }
    }

    list(entropy_rate1 = rate)
  }, error = function(e) default_result)
}

#' Goodness-of-fit to normality
#'
#' Anderson-Darling and Cramér-von Mises statistics.
#'
#' @param x Numeric vector
#' @return Named list with ad_stat, cvm_stat
#' @export
ts_goodness_normality <- function(x) {
  n <- length(x)

  default_result <- list(
    ad_stat = NA_real_,
    cvm_stat = NA_real_
  )

  if (n < 10) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 10) return(default_result)

  tryCatch({
    z <- (x - mean(x)) / sd(x)
    z <- sort(z)
    n <- length(z)

    F_z <- pnorm(z)

    # Anderson-Darling
    i <- 1:n
    ad <- -n - mean((2 * i - 1) * (log(F_z) + log(1 - rev(F_z))))

    # Cramér-von Mises
    cvm <- sum((F_z - (2 * i - 1) / (2 * n))^2) + 1 / (12 * n)

    list(
      ad_stat = ad,
      cvm_stat = cvm
    )
  }, error = function(e) default_result)
}

#' Page-Hinkley test
#'
#' Change point detection statistic and index.
#'
#' @param x Numeric vector
#' @return Named list with ph_stat, ph_cp_index
#' @export
ts_page_hinkley <- function(x) {
  n <- length(x)

  default_result <- list(
    ph_stat = NA_real_,
    ph_cp_index = NA_integer_
  )

  if (n < 10) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 10) return(default_result)

  tryCatch({
    delta <- 0.01 * sd(x)
    mu_hat <- mean(x)

    PH <- cumsum(x - mu_hat - delta)
    m <- cummin(PH)
    diff <- PH - m

    list(
      ph_stat = max(diff),
      ph_cp_index = which.max(diff)
    )
  }, error = function(e) default_result)
}

#' Integrated autocorrelation time
#'
#' Sum of ACF until first nonpositive or threshold.
#'
#' @param x Numeric vector
#' @return Named list with int_autocorr_time
#' @export
ts_iat <- function(x) {
  n <- length(x)

  default_result <- list(int_autocorr_time = NA_real_)

  if (n < 20) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 20) return(default_result)

  tryCatch({
    max_lag <- min(50, floor(n / 2))
    ac <- acf(x, lag.max = max_lag, plot = FALSE)$acf[-1]

    # Sum until first nonpositive or |ρ| < 1/e
    threshold <- 1 / exp(1)
    tau <- 1

    for (k in seq_along(ac)) {
      if (ac[k] <= 0 || abs(ac[k]) < threshold) break
      tau <- tau + 2 * ac[k]
    }

    list(int_autocorr_time = tau)
  }, error = function(e) default_result)
}

#' Peak symmetry statistics
#'
#' Mean and SD of peak symmetry (left vs right area).
#'
#' @param x Numeric vector
#' @return Named list with peak_symmetry_mean, peak_symmetry_sd
#' @export
ts_peak_symmetry <- function(x) {
  n <- length(x)

  default_result <- list(
    peak_symmetry_mean = NA_real_,
    peak_symmetry_sd = NA_real_
  )

  if (n < 20) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 20) return(default_result)

  tryCatch({
    # Find peaks (local maxima)
    peaks <- which(diff(sign(diff(x))) == -2) + 1
    if (length(peaks) < 2) return(default_result)

    # Filter peaks that are too close to edges
    window <- 5
    peaks <- peaks[peaks > window & peaks < (n - window)]
    if (length(peaks) < 2) return(default_result)

    symmetries <- numeric(length(peaks))

    for (i in seq_along(peaks)) {
      pk <- peaks[i]
      left <- sum(x[(pk - window):(pk - 1)])
      right <- sum(x[(pk + 1):(pk + window)])
      total <- left + right
      if (total > 0) {
        symmetries[i] <- (right - left) / total
      } else {
        symmetries[i] <- 0
      }
    }

    list(
      peak_symmetry_mean = mean(symmetries),
      peak_symmetry_sd = sd(symmetries)
    )
  }, error = function(e) default_result)
}

#' Zero crossing rate of differences
#'
#' Rate of sign changes in first differences.
#'
#' @param x Numeric vector
#' @return Named list with zcr_diff
#' @export
ts_zcr_diff <- function(x) {
  n <- length(x)

  default_result <- list(zcr_diff = NA_real_)

  if (n < 3) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 3) return(default_result)

  tryCatch({
    dx <- diff(x)
    zcr <- mean(dx[-1] * dx[-length(dx)] < 0)

    list(zcr_diff = zcr)
  }, error = function(e) default_result)
}

#' Harmonicity index
#'
#' Ratio of top 3 peaks to total power.
#'
#' @param x Numeric vector
#' @return Named list with harmonicity_top3
#' @export
ts_harmonicity_index <- function(x) {
  n <- length(x)

  default_result <- list(harmonicity_top3 = NA_real_)

  if (n < 16) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 16) return(default_result)

  tryCatch({
    sp <- spec.pgram(x, taper = 0.1, detrend = TRUE, plot = FALSE)
    P <- sp$spec

    # Exclude DC component
    P <- P[-1]
    if (length(P) < 3) return(default_result)

    # Get top 3 peaks
    top3_sum <- sum(sort(P, decreasing = TRUE)[1:3])
    total <- sum(P)

    harm <- if (total > 0) top3_sum / total else NA_real_

    list(harmonicity_top3 = harm)
  }, error = function(e) default_result)
}

#' Spectral slope R²
#'
#' R² of log-log spectral slope fit.
#'
#' @param x Numeric vector
#' @return Named list with spectral_slope_r2
#' @export
ts_spectral_slope_r2 <- function(x) {
  n <- length(x)

  default_result <- list(spectral_slope_r2 = NA_real_)

  if (n < 16) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 16) return(default_result)

  tryCatch({
    sp <- spec.pgram(x, taper = 0.1, detrend = TRUE, plot = FALSE)
    P <- sp$spec[sp$spec > 0]
    F <- sp$freq[sp$spec > 0]

    if (length(P) < 5) return(default_result)

    fit <- lm(log(P) ~ F)
    r2 <- summary(fit)$r.squared

    list(spectral_slope_r2 = r2)
  }, error = function(e) default_result)
}

#' Quantile spread trend
#'
#' Slope of rolling IQR over time.
#'
#' @param x Numeric vector
#' @return Named list with qsp_slope
#' @export
ts_quantile_spread_trend <- function(x) {
  n <- length(x)

  default_result <- list(qsp_slope = NA_real_)

  if (n < 30) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 30) return(default_result)

  tryCatch({
    window_size <- max(10, floor(n / 10))
    n_windows <- floor(n / window_size)

    if (n_windows < 3) return(default_result)

    iqrs <- numeric(n_windows)

    for (i in 1:n_windows) {
      start <- (i - 1) * window_size + 1
      end <- min(i * window_size, n)
      window <- x[start:end]
      iqrs[i] <- IQR(window, type = 7)
    }

    fit <- lm(iqrs ~ seq_along(iqrs))
    slope <- coef(fit)[2]

    list(qsp_slope = unname(slope))
  }, error = function(e) default_result)
}

#' Crossings around mean and median
#'
#' Counts of crossings around mean and median.
#'
#' @param x Numeric vector
#' @return Named list with crossings_mean, crossings_median
#' @export
ts_crossings_mean_median <- function(x) {
  n <- length(x)

  default_result <- list(
    crossings_mean = NA_real_,
    crossings_median = NA_real_
  )

  if (n < 3) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 3) return(default_result)

  tryCatch({
    m <- mean(x)
    med <- median(x)

    cross_mean <- sum((x[-n] - m) * (x[-1] - m) < 0) / (n - 1)
    cross_median <- sum((x[-n] - med) * (x[-1] - med) < 0) / (n - 1)

    list(
      crossings_mean = cross_mean,
      crossings_median = cross_median
    )
  }, error = function(e) default_result)
}

#' Normalized zero-lag autocovariance
#'
#' Ratio of variance to mean squared.
#'
#' @param x Numeric vector
#' @return Named list with acov0_norm
#' @export
ts_autocov_zero_lag_norm <- function(x) {
  n <- length(x)

  default_result <- list(acov0_norm = NA_real_)

  if (n < 2) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 2) return(default_result)

  tryCatch({
    v <- var(x)
    m2 <- mean(x^2)

    acov0 <- if (m2 > 0) v / m2 else NA_real_

    list(acov0_norm = acov0)
  }, error = function(e) default_result)
}

#' Event rate over threshold
#'
#' Rate and mean exceedance over threshold.
#'
#' @param x Numeric vector
#' @return Named list with event_rate_t, mean_exceedance_t
#' @export
ts_event_rate_over_threshold <- function(x) {
  n <- length(x)

  default_result <- list(
    event_rate_t = NA_real_,
    mean_exceedance_t = NA_real_
  )

  if (n < 10) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 10) return(default_result)

  tryCatch({
    thr <- median(x) + sd(x)
    events <- x > thr
    event_rate <- mean(events)

    exceedances <- x[events] - thr
    mean_exc <- if (length(exceedances) > 0) mean(exceedances) else 0

    list(
      event_rate_t = event_rate,
      mean_exceedance_t = mean_exc
    )
  }, error = function(e) default_result)
}

#' Geweke-Porter-Hudak long-memory estimate
#'
#' Estimates long-range dependence parameter via periodogram regression.
#'
#' @param x Numeric vector
#' @return Named list with gph_d
#' @export
ts_gph_d <- function(x) {
  n <- length(x)

  default_result <- list(gph_d = NA_real_)

  if (n < 50) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 50) return(default_result)

  tryCatch({
    # Compute periodogram
    sp <- spec.pgram(x - mean(x), taper = 0, detrend = FALSE, plot = FALSE, fast = FALSE)
    P <- sp$spec
    freq <- sp$freq

    # Use low frequencies (typically m = n^0.5 to n^0.6)
    m <- floor(length(P)^0.5)
    m <- max(5, min(m, floor(length(P) / 2)))

    # Log-log regression near zero frequency
    log_P <- log(P[1:m])
    log_freq <- log(freq[1:m])

    if (any(!is.finite(log_P)) || any(!is.finite(log_freq))) return(default_result)

    fit <- lm(log_P ~ log_freq)
    d_est <- -(coef(fit)[2] + 1) / 2

    list(gph_d = unname(d_est))
  }, error = function(e) default_result)
}

#' Mean Absolute Scaled Error
#'
#' In-sample forecast error normalized by naive forecast MAE.
#'
#' @param x Numeric vector
#' @return Named list with mase, smase
#' @export
ts_mase <- function(x) {
  n <- length(x)

  default_result <- list(
    mase = NA_real_,
    smase = NA_real_
  )

  if (n < 20) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 20) return(default_result)

  tryCatch({
    # Simple in-sample one-step-ahead forecast using mean
    forecast_vals <- rep(mean(x[-n]), n - 1)
    errors <- x[-1] - forecast_vals

    # Naive forecast baseline (random walk)
    naive_errors <- diff(x)
    mae_naive <- mean(abs(naive_errors))

    # MASE
    mase <- if (mae_naive > 0) mean(abs(errors)) / mae_naive else NA_real_

    # Seasonal MASE (using detected period)
    freq_result <- ts_frequency(x)
    period <- as.integer(freq_result$freq_est)

    if (period > 1 && length(x) > 2 * period) {
      seasonal_naive_errors <- x[(period + 1):n] - x[1:(n - period)]
      mae_seasonal_naive <- mean(abs(seasonal_naive_errors))
      smase <- if (mae_seasonal_naive > 0) mean(abs(errors)) / mae_seasonal_naive else NA_real_
    } else {
      smase <- NA_real_
    }

    list(
      mase = mase,
      smase = smase
    )
  }, error = function(e) default_result)
}

#' BDS test for nonlinearity
#'
#' Correlation integral-based test for independence.
#'
#' @param x Numeric vector
#' @return Named list with bds_stat, bds_pvalue
#' @export
ts_bds <- function(x) {
  n <- length(x)

  default_result <- list(
    bds_stat = NA_real_,
    bds_pvalue = NA_real_
  )

  if (n < 50) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 50) return(default_result)

  tryCatch({
    # Standardize
    z <- (x - mean(x)) / sd(x)

    # Embedding dimension and threshold
    m <- 2
    epsilon <- 0.5 * sd(z)

    # Create embedded vectors
    n_embed <- n - m + 1
    embedded <- matrix(0, n_embed, m)
    for (i in 1:m) {
      embedded[, i] <- z[i:(n - m + i)]
    }

    # Compute correlation integrals
    count_m <- 0
    count_1 <- 0

    for (i in 1:(n_embed - 1)) {
      for (j in (i + 1):n_embed) {
        # m-dimensional distance
        dist_m <- max(abs(embedded[i, ] - embedded[j, ]))
        if (dist_m < epsilon) count_m <- count_m + 1

        # 1-dimensional distance
        dist_1 <- abs(z[i] - z[j])
        if (dist_1 < epsilon) count_1 <- count_1 + 1
      }
    }

    n_pairs <- n_embed * (n_embed - 1) / 2

    C_m <- count_m / n_pairs
    C_1 <- count_1 / n_pairs

    # BDS statistic
    V <- sqrt(n_embed) * (C_m - C_1^m)

    # Approximate variance (simplified)
    var_approx <- m * (C_1^(m - 1))^2 * C_1 * (1 - C_1)

    if (var_approx > 0) {
      bds_stat <- V / sqrt(var_approx)
      bds_pval <- 2 * (1 - pnorm(abs(bds_stat)))
    } else {
      bds_stat <- NA_real_
      bds_pval <- NA_real_
    }

    list(
      bds_stat = bds_stat,
      bds_pvalue = bds_pval
    )
  }, error = function(e) default_result)
}

#' Dispersion entropy
#'
#' Robust entropy using value-to-class mapping.
#'
#' @param x Numeric vector
#' @return Named list with dispersion_entropy
#' @export
ts_entropy_dispersion <- function(x) {
  n <- length(x)

  default_result <- list(dispersion_entropy = NA_real_)

  if (n < 20) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 20) return(default_result)

  tryCatch({
    # Parameters
    m <- 3  # Embedding dimension
    c <- 5  # Number of classes

    # Map to classes using normal CDF
    z <- (x - mean(x)) / sd(x)
    mapped <- pnorm(z)
    classes <- ceiling(mapped * c)
    classes[classes < 1] <- 1
    classes[classes > c] <- c

    # Create dispersion patterns
    n_patterns <- n - m + 1
    patterns <- character(n_patterns)

    for (i in 1:n_patterns) {
      patterns[i] <- paste(classes[i:(i + m - 1)], collapse = "-")
    }

    # Count pattern frequencies
    pattern_counts <- table(patterns)
    probs <- pattern_counts / n_patterns

    # Compute entropy
    disp_ent <- -sum(probs * log(probs))

    # Normalize by max entropy
    max_ent <- log(c^m)
    disp_ent_norm <- if (max_ent > 0) disp_ent / max_ent else NA_real_

    list(dispersion_entropy = disp_ent_norm)
  }, error = function(e) default_result)
}

#' Multiscale permutation entropy
#'
#' Ordinal-pattern entropy over coarse-grained scales.
#'
#' @param x Numeric vector
#' @return Named list with mpe_2, mpe_3, mpe_5
#' @export
ts_multiscale_perm <- function(x) {
  n <- length(x)

  default_result <- list(
    mpe_2 = NA_real_,
    mpe_3 = NA_real_,
    mpe_5 = NA_real_
  )

  if (n < 50) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 50) return(default_result)

  # Helper function to compute permutation entropy
  perm_entropy <- function(y, m = 3) {
    n_y <- length(y)
    if (n_y < m + 1) return(NA_real_)

    n_patterns <- n_y - m + 1
    patterns <- character(n_patterns)

    for (i in 1:n_patterns) {
      segment <- y[i:(i + m - 1)]
      pattern <- paste(order(segment), collapse = "-")
      patterns[i] <- pattern
    }

    counts <- table(patterns)
    probs <- counts / n_patterns

    ent <- -sum(probs * log(probs))

    # Normalize by max entropy (log(m!))
    max_ent <- log(factorial(m))
    if (max_ent > 0) ent / max_ent else NA_real_
  }

  # Coarse-graining function
  coarse_grain <- function(y, scale) {
    n_y <- length(y)
    n_cg <- floor(n_y / scale)
    cg <- numeric(n_cg)

    for (i in 1:n_cg) {
      start <- (i - 1) * scale + 1
      end <- i * scale
      cg[i] <- mean(y[start:end])
    }

    cg
  }

  tryCatch({
    # Scale 2
    cg2 <- coarse_grain(x, 2)
    mpe2 <- perm_entropy(cg2, m = 3)

    # Scale 3
    cg3 <- coarse_grain(x, 3)
    mpe3 <- perm_entropy(cg3, m = 3)

    # Scale 5
    cg5 <- coarse_grain(x, 5)
    mpe5 <- perm_entropy(cg5, m = 3)

    list(
      mpe_2 = mpe2,
      mpe_3 = mpe3,
      mpe_5 = mpe5
    )
  }, error = function(e) default_result)
}

#' First harmonic fit quality
#'
#' R² and RMSE of first harmonic fit at detected period.
#'
#' @param x Numeric vector
#' @return Named list with harmonic_r2, harmonic_rmse
#' @export
ts_fourier1_resid <- function(x) {
  n <- length(x)

  default_result <- list(
    harmonic_r2 = NA_real_,
    harmonic_rmse = NA_real_
  )

  if (n < 30) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 30) return(default_result)

  tryCatch({
    # Get period estimate
    freq_result <- ts_frequency(x)
    P <- max(2L, as.integer(freq_result$freq_est))

    if (P < 2 || n < 2 * P) return(default_result)

    # Fit first harmonic
    t <- seq_len(n)
    fit <- lm(x ~ I(sin(2 * pi * t / P)) + I(cos(2 * pi * t / P)))

    # Extract R² and RMSE
    r2 <- summary(fit)$r.squared
    residuals <- residuals(fit)
    rmse <- sqrt(mean(residuals^2))

    list(
      harmonic_r2 = r2,
      harmonic_rmse = rmse
    )
  }, error = function(e) default_result)
}

#' Trimmed statistics
#'
#' Trimmed mean and SD at 5% level.
#'
#' @param x Numeric vector
#' @return Named list with trimmed_mean_p05, trimmed_sd_p05
#' @export
ts_trimmed_stats <- function(x) {
  n <- length(x)

  default_result <- list(
    trimmed_mean_p05 = NA_real_,
    trimmed_sd_p05 = NA_real_
  )

  if (n < 10) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 10) return(default_result)

  tryCatch({
    # 5% trimmed mean
    trim_mean <- mean(x, trim = 0.05)

    # 5% trimmed SD
    trim_n <- floor(length(x) * 0.05)
    sorted_x <- sort(x)
    trimmed_x <- sorted_x[(trim_n + 1):(length(x) - trim_n)]
    trim_sd <- sd(trimmed_x)

    list(
      trimmed_mean_p05 = trim_mean,
      trimmed_sd_p05 = trim_sd
    )
  }, error = function(e) default_result)
}

#' Tail asymmetry
#'
#' Upper and lower tail ratios.
#'
#' @param x Numeric vector
#' @return Named list with upper_tail_ratio, lower_tail_ratio
#' @export
ts_tail_asymmetry <- function(x) {
  n <- length(x)

  default_result <- list(
    upper_tail_ratio = NA_real_,
    lower_tail_ratio = NA_real_
  )

  if (n < 20) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 20) return(default_result)

  tryCatch({
    q <- quantile(x, c(0.01, 0.05, 0.95, 0.99), type = 7)

    eps <- 1e-10

    # Upper tail: (Q99 - Q95) / (abs(Q95) + eps)
    upper_ratio <- (q[4] - q[3]) / (abs(q[3]) + eps)

    # Lower tail: (Q05 - Q01) / (abs(Q05) + eps)
    lower_ratio <- (q[2] - q[1]) / (abs(q[2]) + eps)

    list(
      upper_tail_ratio = upper_ratio,
      lower_tail_ratio = lower_ratio
    )
  }, error = function(e) default_result)
}

#' Multiple peaks in PSD
#'
#' Count and contrast of spectral peaks.
#'
#' @param x Numeric vector
#' @return Named list with psd_num_peaks, psd_peak_contrast
#' @export
ts_multipeak_psd <- function(x) {
  n <- length(x)

  default_result <- list(
    psd_num_peaks = NA_real_,
    psd_peak_contrast = NA_real_
  )

  if (n < 16) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 16) return(default_result)

  tryCatch({
    sp <- spec.pgram(x, taper = 0.1, detrend = TRUE, plot = FALSE)
    P <- sp$spec

    # Find local maxima
    n_p <- length(P)
    if (n_p < 3) return(default_result)

    is_peak <- logical(n_p)
    for (i in 2:(n_p - 1)) {
      if (P[i] > P[i - 1] && P[i] > P[i + 1]) {
        is_peak[i] <- TRUE
      }
    }

    # Count peaks above median * factor
    factor <- 1.5
    threshold <- median(P) * factor
    num_peaks <- sum(P[is_peak] > threshold)

    # Contrast: max(P) / median(P)
    contrast <- max(P) / (median(P) + 1e-10)

    list(
      psd_num_peaks = num_peaks,
      psd_peak_contrast = contrast
    )
  }, error = function(e) default_result)
}

#' Block variance ratio
#'
#' Ratio of max to min variance across blocks.
#'
#' @param x Numeric vector
#' @return Named list with var_ratio_blocks
#' @export
ts_block_var_ratio <- function(x) {
  n <- length(x)

  default_result <- list(var_ratio_blocks = NA_real_)

  if (n < 30) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 30) return(default_result)

  tryCatch({
    # Split into k blocks
    k <- max(3, min(10, floor(n / 10)))
    block_size <- floor(n / k)

    block_vars <- numeric(k)

    for (i in 1:k) {
      start <- (i - 1) * block_size + 1
      end <- if (i == k) n else i * block_size
      block <- x[start:end]
      block_vars[i] <- var(block)
    }

    # Filter out zero variances
    block_vars <- block_vars[block_vars > 0]

    if (length(block_vars) < 2) return(default_result)

    var_ratio <- max(block_vars) / min(block_vars)

    list(var_ratio_blocks = var_ratio)
  }, error = function(e) default_result)
}

#' Difference statistics
#'
#' Descriptors on first differences.
#'
#' @param x Numeric vector
#' @return Named list with mean_abs_diff, var_diff, kurt_diff
#' @export
ts_diff_stats <- function(x) {
  n <- length(x)

  default_result <- list(
    mean_abs_diff = NA_real_,
    var_diff = NA_real_,
    kurt_diff = NA_real_
  )

  if (n < 5) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 5) return(default_result)

  tryCatch({
    dx <- diff(x)

    mad <- mean(abs(dx))
    vard <- var(dx)
    kurtd <- if (length(dx) >= 4) cpp_kurtosis(dx) else NA_real_

    list(
      mean_abs_diff = mad,
      var_diff = vard,
      kurt_diff = kurtd
    )
  }, error = function(e) default_result)
}

#' Return statistics
#'
#' Descriptors on returns (relative changes).
#'
#' @param x Numeric vector
#' @return Named list with mean_abs_ret, vol_ret, skew_ret
#' @export
ts_return_stats <- function(x) {
  n <- length(x)

  default_result <- list(
    mean_abs_ret = NA_real_,
    vol_ret = NA_real_,
    skew_ret = NA_real_
  )

  if (n < 5) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 5) return(default_result)

  tryCatch({
    # Compute returns: diff(x) / lag(x), guarding against zeros
    x_lag <- x[-n]
    x_lead <- x[-1]

    # Guard against zeros
    eps <- 1e-10
    returns <- (x_lead - x_lag) / (abs(x_lag) + eps)

    # Remove extreme values (likely from near-zero denominators)
    returns <- returns[abs(returns) < 100]

    if (length(returns) < 3) return(default_result)

    mar <- mean(abs(returns))
    vol <- sd(returns)
    skew <- if (length(returns) >= 3) cpp_skewness(returns) else NA_real_

    list(
      mean_abs_ret = mar,
      vol_ret = vol,
      skew_ret = skew
    )
  }, error = function(e) default_result)
}

#' NA gap statistics
#'
#' NA rate and longest NA run (computed on original data before NA removal).
#'
#' @param x Numeric vector (with potential NAs)
#' @return Named list with na_rate, longest_na_run
#' @export
ts_na_gaps <- function(x) {
  n <- length(x)

  default_result <- list(
    na_rate = NA_real_,
    longest_na_run = NA_integer_
  )

  if (n < 1) return(default_result)

  tryCatch({
    na_mask <- is.na(x)
    na_rate <- mean(na_mask)

    # Find longest NA run using rle
    runs <- rle(na_mask)
    na_runs <- runs$lengths[runs$values]

    longest_run <- if (length(na_runs) > 0) max(na_runs) else 0L

    list(
      na_rate = na_rate,
      longest_na_run = as.integer(longest_run)
    )
  }, error = function(e) default_result)
}

#' Constant segments
#'
#' Count and max length of constant segments.
#'
#' @param x Numeric vector
#' @return Named list with const_seg_count, const_seg_max_len
#' @export
ts_const_segments <- function(x) {
  n <- length(x)

  default_result <- list(
    const_seg_count = NA_integer_,
    const_seg_max_len = NA_integer_
  )

  if (n < 2) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 2) return(default_result)

  tryCatch({
    # Round to avoid floating point issues
    eps <- 1e-10
    x_rounded <- round(x / eps) * eps

    # Identify where diff == 0
    dx <- diff(x_rounded)
    is_const <- abs(dx) < eps * 10

    # Use rle to find runs
    runs <- rle(is_const)
    const_runs <- runs$lengths[runs$values]

    count <- length(const_runs)
    max_len <- if (count > 0) max(const_runs) + 1L else 0L  # +1 because diff reduces length by 1

    list(
      const_seg_count = as.integer(count),
      const_seg_max_len = as.integer(max_len)
    )
  }, error = function(e) default_result)
}

#' Value range over standard deviation
#'
#' Range normalized by standard deviation.
#'
#' @param x Numeric vector
#' @return Named list with range_over_std
#' @export
ts_value_range <- function(x) {
  n <- length(x)

  default_result <- list(range_over_std = NA_real_)

  if (n < 2) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 2) return(default_result)

  tryCatch({
    r <- max(x) - min(x)
    s <- sd(x)

    eps <- 1e-10
    ratio <- r / (s + eps)

    list(range_over_std = ratio)
  }, error = function(e) default_result)
}

#' Pettitt change point test
#'
#' Non-parametric test for change point in distribution.
#'
#' @param x Numeric vector
#' @return Named list with pettitt_stat, pettitt_pval, pettitt_cp_index
#' @export
ts_pettitt <- function(x) {
  n <- length(x)

  default_result <- list(
    pettitt_stat = NA_real_,
    pettitt_pval = NA_real_,
    pettitt_cp_index = NA_integer_
  )

  if (n < 20) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 20) return(default_result)

  tryCatch({
    # Compute U statistic at each possible split point
    U <- numeric(n)

    for (t in 1:(n - 1)) {
      # Mann-Whitney U-like statistic
      n1 <- t
      n2 <- n - t

      x1 <- x[1:t]
      x2 <- x[(t + 1):n]

      # Count how many pairs (i,j) with i in x1, j in x2 have x1[i] > x2[j]
      U_t <- 0
      for (i in 1:n1) {
        for (j in 1:n2) {
          U_t <- U_t + sign(x1[i] - x2[j])
        }
      }

      U[t] <- abs(U_t)
    }

    # Pettitt statistic is max |U|
    K <- max(U)
    cp_idx <- which.max(U)

    # Approximate p-value (valid for large n)
    p_val <- 2 * exp(-6 * K^2 / (n^3 + n^2))

    list(
      pettitt_stat = K,
      pettitt_pval = p_val,
      pettitt_cp_index = as.integer(cp_idx)
    )
  }, error = function(e) default_result)
}

#' Rolling window stability measures
#'
#' Slopes of rolling IQR, variance, and ACF(1).
#'
#' @param x Numeric vector
#' @return Named list with rolling_iqr_slope, rolling_var_slope, rolling_acf1_slope, rolling_cv
#' @export
ts_rolling_stability <- function(x) {
  n <- length(x)

  default_result <- list(
    rolling_iqr_slope = NA_real_,
    rolling_var_slope = NA_real_,
    rolling_acf1_slope = NA_real_,
    rolling_cv = NA_real_
  )

  if (n < 50) return(default_result)

  x <- x[!is.na(x)]
  if (length(x) < 50) return(default_result)

  tryCatch({
    # Window size
    w <- max(10, floor(n / 10))
    n_windows <- floor(n / w)

    if (n_windows < 5) return(default_result)

    iqr_vals <- numeric(n_windows)
    var_vals <- numeric(n_windows)
    acf1_vals <- numeric(n_windows)

    for (i in 1:n_windows) {
      start <- (i - 1) * w + 1
      end <- min(i * w, n)
      window <- x[start:end]

      iqr_vals[i] <- IQR(window, type = 7)
      var_vals[i] <- var(window)

      # ACF(1)
      if (length(window) > 5) {
        ac <- acf(window, lag.max = 1, plot = FALSE)$acf
        acf1_vals[i] <- if (length(ac) > 1) ac[2] else NA_real_
      } else {
        acf1_vals[i] <- NA_real_
      }
    }

    # Fit linear trends
    time_idx <- seq_along(iqr_vals)

    fit_iqr <- lm(iqr_vals ~ time_idx)
    iqr_slope <- coef(fit_iqr)[2]

    fit_var <- lm(var_vals ~ time_idx)
    var_slope <- coef(fit_var)[2]

    valid_acf <- !is.na(acf1_vals)
    if (sum(valid_acf) >= 3) {
      fit_acf <- lm(acf1_vals[valid_acf] ~ time_idx[valid_acf])
      acf1_slope <- coef(fit_acf)[2]
    } else {
      acf1_slope <- NA_real_
    }

    # Coefficient of variation of rolling variances
    cv <- sd(var_vals) / (mean(var_vals) + 1e-10)

    list(
      rolling_iqr_slope = unname(iqr_slope),
      rolling_var_slope = unname(var_slope),
      rolling_acf1_slope = unname(acf1_slope),
      rolling_cv = cv
    )
  }, error = function(e) default_result)
}

#' Cross-correlation function summary
#'
#' Maximum absolute CCF and its lag (bivariate feature).
#'
#' @param x Numeric vector (first series)
#' @param y Numeric vector (second series)
#' @return Named list with ccf_max_abs, ccf_argmax
#' @export
ts_ccf_summary <- function(x, y) {
  n_x <- length(x)
  n_y <- length(y)

  default_result <- list(
    ccf_max_abs = NA_real_,
    ccf_argmax = NA_integer_
  )

  if (n_x < 20 || n_y < 20) return(default_result)

  x <- x[!is.na(x)]
  y <- y[!is.na(y)]

  if (length(x) < 20 || length(y) < 20) return(default_result)
  if (length(x) != length(y)) {
    warning("x and y must have same length after NA removal")
    return(default_result)
  }

  tryCatch({
    # Compute CCF
    max_lag <- min(20, floor(length(x) / 4))
    ccf_result <- ccf(x, y, lag.max = max_lag, plot = FALSE)

    ccf_vals <- as.numeric(ccf_result$acf)
    lags <- ccf_result$lag

    # Find max absolute CCF
    max_abs_idx <- which.max(abs(ccf_vals))
    max_abs_ccf <- abs(ccf_vals[max_abs_idx])
    lag_at_max <- lags[max_abs_idx]

    list(
      ccf_max_abs = max_abs_ccf,
      ccf_argmax = as.integer(lag_at_max)
    )
  }, error = function(e) default_result)
}

#' Simple coherence analysis
#'
#' Coherence between two series (bivariate feature).
#'
#' @param x Numeric vector (first series)
#' @param y Numeric vector (second series)
#' @return Named list with coh_peak, coh_mean
#' @export
ts_coherence_simple <- function(x, y) {
  n_x <- length(x)
  n_y <- length(y)

  default_result <- list(
    coh_peak = NA_real_,
    coh_mean = NA_real_
  )

  if (n_x < 30 || n_y < 30) return(default_result)

  x <- x[!is.na(x)]
  y <- y[!is.na(y)]

  if (length(x) < 30 || length(y) < 30) return(default_result)
  if (length(x) != length(y)) {
    warning("x and y must have same length after NA removal")
    return(default_result)
  }

  tryCatch({
    # Demean
    x <- x - mean(x)
    y <- y - mean(y)

    # Compute spectra
    sp_x <- spec.pgram(x, taper = 0.1, detrend = FALSE, plot = FALSE)
    sp_y <- spec.pgram(y, taper = 0.1, detrend = FALSE, plot = FALSE)

    # Cross-spectrum (simplified: use FFT)
    n <- length(x)
    fft_x <- fft(x)
    fft_y <- fft(y)

    # Cross-spectrum
    cross_spec <- fft_x * Conj(fft_y)

    # Take only positive frequencies
    n_freq <- floor(n / 2)
    cross_spec <- cross_spec[1:n_freq]
    Pxx <- sp_x$spec[1:n_freq]
    Pyy <- sp_y$spec[1:n_freq]

    # Coherence: |Sxy|^2 / (Sxx * Syy)
    coherence <- abs(cross_spec)^2 / (Pxx * Pyy + 1e-10)
    coherence <- pmin(coherence, 1)  # Cap at 1

    list(
      coh_peak = max(coherence, na.rm = TRUE),
      coh_mean = mean(coherence, na.rm = TRUE)
    )
  }, error = function(e) default_result)
}

#' Granger causality (lite version)
#'
#' AIC comparison of AR vs ARX models (bivariate feature).
#'
#' @param x Numeric vector (predictor series)
#' @param y Numeric vector (response series)
#' @return Named list with granger_aic_diff, granger_pval
#' @export
ts_granger_lite <- function(x, y) {
  n_x <- length(x)
  n_y <- length(y)

  default_result <- list(
    granger_aic_diff = NA_real_,
    granger_pval = NA_real_
  )

  if (n_x < 30 || n_y < 30) return(default_result)

  x <- x[!is.na(x)]
  y <- y[!is.na(y)]

  if (length(x) < 30 || length(y) < 30) return(default_result)
  if (length(x) != length(y)) {
    warning("x and y must have same length after NA removal")
    return(default_result)
  }

  tryCatch({
    n <- length(y)
    p <- 2  # Lag order

    # Create lagged matrices
    y_lag <- embed(y, p + 1)
    x_lag <- embed(x, p + 1)[, -1]  # Remove contemporaneous x

    # Align lengths
    n_eff <- nrow(y_lag)
    y_current <- y_lag[, 1]
    y_lagged <- y_lag[, -1]

    # AR model: y ~ y_{t-1}, y_{t-2}
    fit_ar <- lm(y_current ~ y_lagged)
    aic_ar <- AIC(fit_ar)

    # ARX model: y ~ y_{t-1}, y_{t-2}, x_{t-1}, x_{t-2}
    fit_arx <- lm(y_current ~ y_lagged + x_lag)
    aic_arx <- AIC(fit_arx)

    # AIC difference (negative means ARX is better)
    aic_diff <- aic_arx - aic_ar

    # F-test for significance
    anova_result <- anova(fit_ar, fit_arx)
    p_val <- anova_result$`Pr(>F)`[2]

    list(
      granger_aic_diff = aic_diff,
      granger_pval = p_val
    )
  }, error = function(e) default_result)
}
