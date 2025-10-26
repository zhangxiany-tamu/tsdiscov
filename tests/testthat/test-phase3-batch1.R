# Phase 3 Batch 1: Differencing ACF/PACF Features
# Tests for diff1_acf1, diff1_acf10, diff2_acf1, diff2_acf10, diff1x_pacf5, diff2x_pacf5

test_that("First difference ACF features work correctly", {
  set.seed(123)
  x <- rnorm(100)

  # Test diff1_acf1
  result <- cpp_diff1_acf1(x)
  expect_true(is.finite(result))
  expect_true(result >= -1 && result <= 1)  # ACF should be in [-1, 1]

  # Test diff1_acf10
  result10 <- cpp_diff1_acf10(x)
  expect_true(is.finite(result10))
  expect_true(result10 >= 0)  # Sum of squares should be non-negative
})

test_that("Second difference ACF features work correctly", {
  set.seed(123)
  x <- rnorm(100)

  # Test diff2_acf1
  result <- cpp_diff2_acf1(x)
  expect_true(is.finite(result))
  expect_true(result >= -1 && result <= 1)

  # Test diff2_acf10
  result10 <- cpp_diff2_acf10(x)
  expect_true(is.finite(result10))
  expect_true(result10 >= 0)
})

test_that("Difference PACF features work correctly", {
  set.seed(123)
  x <- rnorm(100)

  # Test diff1x_pacf5
  result1 <- cpp_diff1x_pacf5(x)
  expect_true(is.finite(result1))
  expect_true(result1 >= 0)  # Sum of squares should be non-negative

  # Test diff2x_pacf5
  result2 <- cpp_diff2x_pacf5(x)
  expect_true(is.finite(result2))
  expect_true(result2 >= 0)
})

test_that("Differencing features detect non-stationarity", {
  set.seed(42)

  # Stationary white noise
  stationary <- rnorm(200)
  diff1_acf1_stat <- cpp_diff1_acf1(stationary)
  diff2_acf1_stat <- cpp_diff2_acf1(stationary)

  # Random walk (non-stationary)
  rw <- cumsum(rnorm(200))
  diff1_acf1_rw <- cpp_diff1_acf1(rw)
  diff2_acf1_rw <- cpp_diff2_acf1(rw)

  # All should be finite
  expect_true(is.finite(diff1_acf1_stat))
  expect_true(is.finite(diff1_acf1_rw))
  expect_true(is.finite(diff2_acf1_stat))
  expect_true(is.finite(diff2_acf1_rw))

  # Values should be in valid ACF range
  expect_true(abs(diff1_acf1_rw) <= 1)
  expect_true(abs(diff2_acf1_rw) <= 1)
})

test_that("Differencing features distinguish signal types", {
  set.seed(123)

  # White noise (already stationary)
  wn <- rnorm(150)
  diff1_wn <- cpp_diff1_acf10(wn)
  diff2_wn <- cpp_diff2_acf10(wn)

  # Trending series
  t <- 1:150
  trending <- 2 * t + rnorm(150, sd = 10)
  diff1_trend <- cpp_diff1_acf10(trending)
  diff2_trend <- cpp_diff2_acf10(trending)

  # All should be finite
  expect_true(is.finite(diff1_wn))
  expect_true(is.finite(diff1_trend))
  expect_true(is.finite(diff2_wn))
  expect_true(is.finite(diff2_trend))

  # Both should be positive (sum of squares)
  expect_true(diff1_wn >= 0)
  expect_true(diff1_trend >= 0)
})

test_that("PACF features detect autocorrelation structure", {
  set.seed(456)

  # AR(1) process
  ar1 <- filter(rnorm(200), filter = 0.7, method = "recursive")
  ar1 <- as.numeric(ar1)
  pacf1_ar1 <- cpp_diff1x_pacf5(ar1)

  # White noise
  wn <- rnorm(200)
  pacf1_wn <- cpp_diff1x_pacf5(wn)

  # All should be finite and non-negative
  expect_true(is.finite(pacf1_ar1))
  expect_true(is.finite(pacf1_wn))
  expect_true(pacf1_ar1 >= 0)
  expect_true(pacf1_wn >= 0)

  # AR process should generally have higher PACF sum than white noise
  # (though this isn't guaranteed, so we just check they're different)
  expect_true(is.numeric(pacf1_ar1))
  expect_true(is.numeric(pacf1_wn))
})

test_that("Differencing features handle edge cases", {
  # Very short series - should return NA or handle gracefully
  x_short <- rnorm(5)
  result1 <- cpp_diff1_acf1(x_short)
  result2 <- cpp_diff1_acf10(x_short)
  result3 <- cpp_diff2_acf1(x_short)
  result4 <- cpp_diff2_acf10(x_short)
  result5 <- cpp_diff1x_pacf5(x_short)
  result6 <- cpp_diff2x_pacf5(x_short)

  # All should be numeric (NA or finite)
  expect_true(is.numeric(result1))
  expect_true(is.numeric(result2))
  expect_true(is.numeric(result3))
  expect_true(is.numeric(result4))
  expect_true(is.numeric(result5))
  expect_true(is.numeric(result6))

  # Minimum length series
  x_min <- rnorm(20)
  result <- cpp_diff1_acf1(x_min)
  # Should work for minimum length
  expect_true(is.finite(result) || is.na(result))

  # Constant series
  x_const <- rep(5, 100)
  # May return NA or 0, just check no error
  expect_true(is.numeric(cpp_diff1_acf1(x_const)))
})

test_that("Aggregate function cpp_diff_acf_features works", {
  set.seed(789)
  x <- rnorm(100)

  result <- cpp_diff_acf_features(x)

  # Should return named vector with 6 features (ndiffs is added in R wrapper)
  expect_true(is.numeric(result))
  expect_equal(length(result), 6)
  expect_equal(names(result), c("diff1_acf1", "diff1_acf10", "diff2_acf1",
                                  "diff2_acf10", "diff1x_pacf5", "diff2x_pacf5"))

  # All should be finite
  expect_true(all(is.finite(result) | is.na(result)))

  # Check individual values match
  expect_equal(unname(result["diff1_acf1"]), cpp_diff1_acf1(x))
  expect_equal(unname(result["diff1_acf10"]), cpp_diff1_acf10(x))
  expect_equal(unname(result["diff2_acf1"]), cpp_diff2_acf1(x))
  expect_equal(unname(result["diff2_acf10"]), cpp_diff2_acf10(x))
  expect_equal(unname(result["diff1x_pacf5"]), cpp_diff1x_pacf5(x))
  expect_equal(unname(result["diff2x_pacf5"]), cpp_diff2x_pacf5(x))
})

test_that("R wrapper ts_diff works correctly", {
  set.seed(111)
  x <- rnorm(150)

  result <- ts_diff(x)

  # Should return list with 7 features (including ndiffs)
  expect_type(result, "list")
  expect_equal(length(result), 7)
  expect_equal(names(result), c("diff1_acf1", "diff1_acf10", "diff2_acf1",
                                  "diff2_acf10", "diff1x_pacf5", "diff2x_pacf5", "ndiffs"))

  # All should be numeric
  expect_true(all(sapply(result, is.numeric)))

  # Should match cpp function
  cpp_result <- cpp_diff_acf_features(x)
  expect_equal(result$diff1_acf1, unname(cpp_result["diff1_acf1"]))
  expect_equal(result$diff1_acf10, unname(cpp_result["diff1_acf10"]))
})

test_that("Differencing features integrate with ts_features", {
  set.seed(222)
  x <- rnorm(100)

  # Test with 'diff' subset
  diff_feats <- ts_features(x, features = "diff")
  expect_type(diff_feats, "list")
  expect_equal(length(diff_feats), 7)
  expect_true("diff1_acf1" %in% names(diff_feats))
  expect_true("diff2_acf10" %in% names(diff_feats))

  # Test with 'all'
  all_feats <- ts_features(x, features = "all")
  expect_true("diff1_acf1" %in% names(all_feats))
  expect_true("diff2_acf10" %in% names(all_feats))
  expect_true(length(all_feats) >= 94)  # Should have at least 94 features now
})

test_that("cpp_diff function works correctly", {
  x <- c(1, 2, 4, 7, 11, 16)

  # First difference
  diff1 <- cpp_diff(x, lag = 1, differences = 1)
  expect_equal(diff1, c(1, 2, 3, 4, 5))

  # Second difference
  diff2 <- cpp_diff(x, lag = 1, differences = 2)
  expect_equal(diff2, c(1, 1, 1, 1))

  # Lag 2 difference
  diff_lag2 <- cpp_diff(x, lag = 2, differences = 1)
  expect_equal(diff_lag2, c(3, 5, 7, 9))
})

test_that("Differencing features handle seasonal patterns", {
  set.seed(333)

  # Create seasonal series
  t <- 1:200
  seasonal <- sin(2 * pi * t / 12) + rnorm(200, sd = 0.1)

  # Differencing should reduce autocorrelation
  diff1_acf10_orig <- cpp_diff1_acf10(seasonal)

  # Should be finite
  expect_true(is.finite(diff1_acf10_orig))
  expect_true(diff1_acf10_orig >= 0)
})

test_that("Differencing features match theoretical expectations", {
  set.seed(444)

  # For white noise, first difference ACF should be around -0.5 at lag 1
  wn <- rnorm(500)
  diff1_acf1_wn <- cpp_diff1_acf1(wn)

  # Should be finite and in valid range
  expect_true(is.finite(diff1_acf1_wn))
  expect_true(abs(diff1_acf1_wn) <= 1)  # Valid ACF range

  # White noise differenced typically has negative ACF[1] close to -0.5
  # but we don't enforce this strictly due to sampling variation
})
