library(testthat)
library(tsdiscov)

# ==============================================================================
# Test Poincaré Features
# ==============================================================================

test_that("Poincaré: white noise has ratio near 1", {
  set.seed(123)
  x <- rnorm(100)

  result <- ts_poincare(x)

  expect_type(result, "list")
  expect_equal(length(result), 3)
  expect_true(all(c("poincare_sd1", "poincare_sd2", "poincare_ratio") %in% names(result)))

  # For white noise, SD1/SD2 ratio should be near 1
  expect_true(result$poincare_ratio > 0.5 && result$poincare_ratio < 2)
})

test_that("Poincaré: random walk has different structure", {
  set.seed(456)
  x <- cumsum(rnorm(100))

  result <- ts_poincare(x)

  expect_true(all(!is.na(unlist(result))))
  expect_gt(result$poincare_sd2, result$poincare_sd1)
})

test_that("Poincaré: short series returns NA", {
  x <- rnorm(5)

  result <- ts_poincare(x)

  expect_true(all(is.na(unlist(result))))
})

# ==============================================================================
# Test Run Statistics
# ==============================================================================

test_that("Run statistics: white noise has moderate runs", {
  set.seed(789)
  x <- rnorm(200)

  result <- ts_runs(x)

  expect_type(result, "list")
  expect_equal(length(result), 4)
  expect_true(all(c("runs_n", "runs_mean_length", "runs_updown_ratio",
                     "runs_switch_rate") %in% names(result)))

  # White noise should have many runs
  expect_gt(result$runs_n, 50)
  # Ratio should be near 1
  expect_true(result$runs_updown_ratio > 0.5 && result$runs_updown_ratio < 2)
})

test_that("Run statistics: trend has fewer runs", {
  x <- seq(1, 100)

  result <- ts_runs(x)

  # Strong trend should have very few runs
  expect_lt(result$runs_n, 10)
  expect_gt(result$runs_mean_length, 10)
})

test_that("Run statistics: short series returns NA", {
  x <- rnorm(5)

  result <- ts_runs(x)

  expect_true(all(is.na(unlist(result))))
})

# ==============================================================================
# Test Shape Factors
# ==============================================================================

test_that("Shape factors: returns correct structure", {
  set.seed(111)
  x <- rnorm(100)

  result <- ts_shape_factors(x)

  expect_type(result, "list")
  expect_equal(length(result), 4)
  expect_true(all(c("shape_crest", "shape_impulse", "shape_shape",
                     "shape_clearance") %in% names(result)))
  expect_true(all(!is.na(unlist(result))))
})

test_that("Shape factors: sine wave has low crest factor", {
  x <- sin(2 * pi * seq(0, 10, length.out = 200))

  result <- ts_shape_factors(x)

  # Sine wave has crest factor sqrt(2) ≈ 1.414
  expect_true(result$shape_crest > 1.2 && result$shape_crest < 1.6)
})

test_that("Shape factors: impulse has high crest factor", {
  x <- c(rep(0.1, 50), 10, rep(0.1, 49))

  result <- ts_shape_factors(x)

  # Impulse should have very high crest and impulse factors
  expect_gt(result$shape_crest, 5)
  expect_gt(result$shape_impulse, 10)
})

test_that("Shape factors: short series returns NA", {
  x <- rnorm(5)

  result <- ts_shape_factors(x)

  expect_true(all(is.na(unlist(result))))
})

# ==============================================================================
# Test Forecastability
# ==============================================================================

test_that("Forecastability: white noise has low forecastability", {
  set.seed(222)
  x <- rnorm(200)

  result <- ts_forecastability(x)

  expect_type(result, "list")
  expect_equal(length(result), 1)
  expect_true("forecastability" %in% names(result))

  # White noise should have low forecastability
  expect_lt(result$forecastability, 0.5)
})

test_that("Forecastability: periodic signal has high forecastability", {
  x <- sin(2 * pi * seq(0, 10, length.out = 200)) + rnorm(200, sd = 0.1)

  result <- ts_forecastability(x)

  # Periodic signal should have higher forecastability
  expect_gt(result$forecastability, 0.3)
})

test_that("Forecastability: short series returns NA", {
  x <- rnorm(15)

  result <- ts_forecastability(x)

  expect_true(is.na(result$forecastability))
})

# ==============================================================================
# Test Frequency Detection
# ==============================================================================

test_that("Frequency detection: returns correct structure", {
  set.seed(333)
  x <- rnorm(100)

  result <- ts_frequency(x)

  expect_type(result, "list")
  expect_equal(length(result), 2)
  expect_true(all(c("freq_est", "nperiods") %in% names(result)))
  expect_true(is.integer(result$freq_est))
  expect_true(is.integer(result$nperiods))
})

test_that("Frequency detection: detects periodic signal", {
  # Strong periodic signal with period 12
  x <- sin(2 * pi * seq(0, 20, length.out = 240) / 12) + rnorm(240, sd = 0.05)

  result <- ts_frequency(x)

  # Should return a valid frequency estimate
  expect_true(is.integer(result$freq_est))
  expect_gte(result$freq_est, 1)
  # Should detect multiple periods
  expect_gte(result$nperiods, 0)
})

test_that("Frequency detection: white noise has low prominence", {
  set.seed(444)
  x <- rnorm(200)

  result <- ts_frequency(x)

  # White noise might detect some period but typically with low nperiods
  expect_true(is.integer(result$freq_est))
  expect_true(is.integer(result$nperiods))
})

test_that("Frequency detection: short series returns default", {
  x <- rnorm(5)

  result <- ts_frequency(x)

  expect_equal(result$freq_est, 1L)
  expect_equal(result$nperiods, 0L)
})

# ==============================================================================
# Test ARIMA Diagnostics
# ==============================================================================

test_that("ARIMA diagnostics: white noise has reasonable p-values", {
  set.seed(444)
  x <- rnorm(100)

  result <- ts_arima_diag(x)

  expect_type(result, "list")
  expect_equal(length(result), 3)
  expect_true(all(c("arima_ljungbox_pval", "arima_normality_pval",
                     "arima_arch_pval") %in% names(result)))

  # White noise should typically have high Ljung-Box p-value
  expect_gt(result$arima_ljungbox_pval, 0.01)
  # All p-values should be valid probabilities
  expect_true(all(sapply(result, function(p) p >= 0 && p <= 1)))
})

test_that("ARIMA diagnostics: returns all numeric values", {
  set.seed(555)
  x <- arima.sim(list(ar = 0.5), n = 100)

  result <- ts_arima_diag(x)

  expect_true(all(sapply(result, is.numeric)))
})

test_that("ARIMA diagnostics: short series returns NA", {
  x <- rnorm(20)

  result <- ts_arima_diag(x)

  expect_true(all(is.na(unlist(result))))
})

# ==============================================================================
# Test Integration with ts_features
# ==============================================================================

test_that("New features integrate into ts_features()", {
  set.seed(666)
  x <- rnorm(100)

  # Test individual extraction
  poincare_only <- ts_features(x, features = "poincare")
  expect_equal(length(poincare_only), 3)

  runs_only <- ts_features(x, features = "runs")
  expect_equal(length(runs_only), 4)

  shape_only <- ts_features(x, features = "shape_factors")
  expect_equal(length(shape_only), 4)

  forecast_only <- ts_features(x, features = "forecastability")
  expect_equal(length(forecast_only), 1)

  freq_only <- ts_features(x, features = "frequency")
  expect_equal(length(freq_only), 2)  # freq_est and nperiods

  arima_only <- ts_features(x, features = "arima_diag")
  expect_equal(length(arima_only), 3)

  # Test in 'all'
  all_feats <- ts_features(x, features = "all")
  expect_true("poincare_sd1" %in% names(all_feats))
  expect_true("runs_n" %in% names(all_feats))
  expect_true("shape_crest" %in% names(all_feats))
  expect_true("forecastability" %in% names(all_feats))
  expect_true("freq_est" %in% names(all_feats))
  expect_true("nperiods" %in% names(all_feats))
  expect_true("arima_ljungbox_pval" %in% names(all_feats))
  expect_true("pettitt_stat" %in% names(all_feats))
  expect_true("rolling_iqr_slope" %in% names(all_feats))

  # Should have 352 total features
  expect_equal(length(all_feats), EXPECTED_ALL_FEATURES)
})
