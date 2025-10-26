#' Feature Sets Registry
#'
#' Internal registry mapping feature set names to their extraction functions.
#' This registry enables discovery, validation, and programmatic access to all
#' available time series feature sets.
#'
#' @keywords internal
#' @noRd
.ts_feature_sets_registry <- function() {
  list(
    stats = ts_stats,
    acf = ts_acf,
    pacf = ts_pacf,
    trend = ts_trend,
    structure = ts_structure,
    fft = ts_fft,
    scaling = ts_scaling,
    stl = ts_stl,
    diff = ts_diff,
    heterogeneity = ts_heterogeneity,
    stattests = ts_stattests,
    shifts = ts_shifts,
    misc = ts_misc,
    holt = ts_holt,
    compengine = ts_compengine,
    catch22 = ts_catch22,
    entropy = ts_entropy,
    complexity = ts_tsfresh,
    aggregations = ts_tsfresh_supp,
    cwt = ts_cwt,
    basic = ts_basic,
    extrema = ts_extrema,
    ar = ts_ar,
    adf = ts_adf,
    welch = ts_welch,
    matrixprofile = ts_matrix_profile,
    friedrich = ts_friedrich,
    langevin = ts_langevin,
    pacf_extended = ts_pacf_extended,
    permutation_entropy = ts_permutation_entropy,
    trend_extended = ts_trend_extended,
    multiscale_entropy = ts_multiscale_entropy,
    robust_trend = ts_robust_trend,
    recurrence = ts_recurrence,
    ljungbox = ts_ljungbox,
    spectral_shape = ts_spectral_shape,
    poincare = ts_poincare,
    runs = ts_runs,
    shape_factors = ts_shape_factors,
    forecastability = ts_forecastability,
    frequency = ts_frequency,
    arima_diag = ts_arima_diag,
    hjorth = ts_hjorth,
    line_length = ts_line_length,
    ssc = ts_ssc,
    turning_points = ts_turning_points,
    monotonicity = ts_monotonicity,
    qcd = ts_qcd,
    iqr_props = ts_iqr_props,
    gini = ts_gini,
    amdf = ts_amdf,
    acf_integrals = ts_acf_integrals,
    spectral_moments = ts_spectral_moments,
    spectral_psr = ts_spectral_psr,
    tkao = ts_tkao,
    time_reversal_multi = ts_time_reversal_multi,
    stl_simple = ts_stl_simple,
    cusum = ts_cusum,
    spectral_qfactor = ts_spectral_qfactor,
    spectral_bands = ts_spectral_bands,
    spectral_crest = ts_spectral_crest,
    tail_index = ts_tail_index,
    robust_shape = ts_robust_shape,
    acf_summary = ts_acf_summary,
    volatility_cluster = ts_volatility_cluster,
    peak_intervals = ts_peak_intervals,
    dwell_times = ts_dwell_times,
    tail_ratios = ts_tail_ratios,
    winsor_stats = ts_winsor_stats,
    renyi2_entropy = ts_renyi2_entropy,
    fourier1 = ts_fourier1,
    ar1_halflife = ts_ar1_halflife,
    cusum_sq = ts_cusum_sq,
    fractal_katz = ts_fractal_katz,
    fractal_petrosian = ts_fractal_petrosian,
    fractal_higuchi = ts_fractal_higuchi,
    entropy_rate_order1 = ts_entropy_rate_order1,
    goodness_normality = ts_goodness_normality,
    page_hinkley = ts_page_hinkley,
    iat = ts_iat,
    peak_symmetry = ts_peak_symmetry,
    zcr_diff = ts_zcr_diff,
    harmonicity_index = ts_harmonicity_index,
    spectral_slope_r2 = ts_spectral_slope_r2,
    quantile_spread_trend = ts_quantile_spread_trend,
    crossings_mean_median = ts_crossings_mean_median,
    autocov_zero_lag_norm = ts_autocov_zero_lag_norm,
    event_rate_over_threshold = ts_event_rate_over_threshold,
    gph_d = ts_gph_d,
    mase = ts_mase,
    bds = ts_bds,
    entropy_dispersion = ts_entropy_dispersion,
    multiscale_perm = ts_multiscale_perm,
    fourier1_resid = ts_fourier1_resid,
    trimmed_stats = ts_trimmed_stats,
    tail_asymmetry = ts_tail_asymmetry,
    multipeak_psd = ts_multipeak_psd,
    block_var_ratio = ts_block_var_ratio,
    diff_stats = ts_diff_stats,
    return_stats = ts_return_stats,
    na_gaps = ts_na_gaps,
    const_segments = ts_const_segments,
    value_range = ts_value_range,
    pettitt = ts_pettitt,
    rolling_stability = ts_rolling_stability
  )
}

#' Get Available Feature Sets
#'
#' Returns a character vector of all available feature set names.
#'
#' @return Character vector of feature set names
#' @export
#' @examples
#' # List all available feature sets
#' feature_sets <- get_feature_sets()
#' head(feature_sets)
get_feature_sets <- function() {
  names(.ts_feature_sets_registry())
}

#' Validate Feature Set Names
#'
#' Checks if requested feature sets are valid. By default, stops with an error
#' if unknown feature sets are found.
#'
#' @param features Character vector of feature set names to validate
#' @param strict Logical; if TRUE (default), stops on unknown features.
#'   If FALSE, warns and returns only valid features.
#' @return Character vector of valid feature set names (invisibly)
#' @export
#' @examples
#' # Valid features pass through
#' validate_feature_sets(c("stats", "acf"))
#'
#' # Unknown features trigger error (strict=TRUE)
#' \dontrun{
#' validate_feature_sets(c("stats", "unknown"))
#' }
#'
#' # With strict=FALSE, warns and returns valid features
#' validate_feature_sets(c("stats", "unknown"), strict = FALSE)
validate_feature_sets <- function(features, strict = TRUE) {
  valid <- get_feature_sets()
  unknown <- setdiff(features, valid)

  if (length(unknown) > 0) {
    msg <- paste("Unknown feature set(s):", paste(unknown, collapse = ", "))
    if (strict) {
      stop(msg, call. = FALSE)
    } else {
      warning(msg, call. = FALSE)
      features <- intersect(features, valid)
    }
  }

  invisible(features)
}

#' List Feature Names for Given Feature Sets
#'
#' Returns the actual feature names that would be extracted by ts_features()
#' for the given feature sets. Useful for dynamic test expectations and
#' understanding what features will be computed.
#'
#' @param x Numeric vector (time series)
#' @param features Character vector of feature set names, or "all"
#' @return Character vector of feature names
#' @export
#' @examples
#' # Get feature names for specific sets
#' x <- rnorm(100)
#' names_stats <- list_ts_feature_names(x, features = "stats")
#' names_acf <- list_ts_feature_names(x, features = c("stats", "acf"))
#'
#' # Get all feature names
#' all_names <- list_ts_feature_names(x, features = "all")
#' length(all_names)
list_ts_feature_names <- function(x, features = "all") {
  # Expand "all" to all feature sets
  if (length(features) == 1 && features == "all") {
    features <- get_feature_sets()
  }

  # Validate
  features <- validate_feature_sets(features, strict = TRUE)

  # Extract features and return names
  result <- ts_features(x, features = features)
  names(result)
}
