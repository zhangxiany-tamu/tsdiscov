library(testthat)
library(tsdiscov)

# ==============================================================================
# Test Ljung-Box
# ==============================================================================

test_that("Ljung-Box: white noise has high p-values", {
  set.seed(123)
  x <- rnorm(100)

  result <- ts_ljungbox(x)

  expect_type(result, "list")
  expect_equal(length(result), 4)
  expect_true(all(c("ljungbox_stat_10", "ljungbox_pval_10",
                     "ljungbox_stat_20", "ljungbox_pval_20") %in% names(result)))

  # White noise should have high p-values (not reject null of no autocorrelation)
  expect_gt(result$ljungbox_pval_10, 0.05)
})

test_that("Ljung-Box: AR(1) series has low p-values", {
  set.seed(456)
  x <- arima.sim(list(ar = 0.8), n = 100)

  result <- ts_ljungbox(x)

  # AR series should reject null of no autocorrelation
  expect_lt(result$ljungbox_pval_10, 0.05)
  expect_lt(result$ljungbox_pval_20, 0.05)
})

test_that("Ljung-Box: short series returns NA", {
  x <- rnorm(20)

  result <- ts_ljungbox(x)

  expect_true(all(is.na(unlist(result))))
})

# ==============================================================================
# Test Spectral Shape
# ==============================================================================

test_that("Spectral shape: white noise has high flatness", {
  set.seed(789)
  x <- rnorm(100)

  result <- ts_spectral_shape(x)

  expect_type(result, "list")
  expect_equal(length(result), 4)
  expect_true(all(c("spectral_flatness", "spectral_slope",
                     "spectral_rolloff_95", "spectral_bandwidth") %in% names(result)))

  # White noise should have relatively high spectral flatness (flatter spectrum)
  expect_gt(result$spectral_flatness, 0.5)
})

test_that("Spectral shape: sine wave has low flatness", {
  # Pure sine wave has concentrated spectrum (low flatness)
  x <- sin(2 * pi * seq(0, 1, length.out = 100) * 5)

  result <- ts_spectral_shape(x)

  # Sine wave should have low spectral flatness (concentrated spectrum)
  expect_lt(result$spectral_flatness, 0.5)
})

test_that("Spectral shape: returns correct types", {
  set.seed(111)
  x <- rnorm(100)

  result <- ts_spectral_shape(x)

  expect_true(all(sapply(result, is.numeric)))
  expect_true(all(!is.na(unlist(result))))
})

test_that("Spectral shape: short series returns NA", {
  x <- rnorm(15)

  result <- ts_spectral_shape(x)

  expect_true(all(is.na(unlist(result))))
})

# ==============================================================================
# Test Integration with ts_features
# ==============================================================================

test_that("New features integrate into ts_features()", {
  set.seed(222)
  x <- rnorm(100)

  # Test individual extraction
  lb_only <- ts_features(x, features = "ljungbox")
  expect_equal(length(lb_only), 4)

  ss_only <- ts_features(x, features = "spectral_shape")
  expect_equal(length(ss_only), 4)

  # Test in 'all'
  all_feats <- ts_features(x, features = "all")
  expect_true("ljungbox_stat_10" %in% names(all_feats))
  expect_true("spectral_flatness" %in% names(all_feats))

  # Should have 352 total features (210 + 4 ljungbox + 4 spectral + 17 base + 32 advanced + 29 GPT)
  expect_equal(length(all_feats), EXPECTED_ALL_FEATURES)
})
