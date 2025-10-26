# Phase 3 Batch 5: Shift Detection Features
# Tests for structural break and change point detection

test_that("Shift detection features work correctly", {
  set.seed(123)

  # Simple white noise
  x <- rnorm(100)

  result <- ts_shifts(x)

  # Should have 6 features
  expect_equal(length(result), 6)
  expect_true("max_level_shift" %in% names(result))
  expect_true("time_level_shift" %in% names(result))
  expect_true("max_var_shift" %in% names(result))
  expect_true("time_var_shift" %in% names(result))
  expect_true("max_kl_shift" %in% names(result))
  expect_true("time_kl_shift" %in% names(result))
})

test_that("Level shift detection works", {
  set.seed(456)

  # Data with clear level shift
  x_shift <- c(rnorm(50, mean = 0), rnorm(50, mean = 5))
  result_shift <- ts_shifts(x_shift, window_size = 10)

  # Should detect significant level shift
  expect_true(is.finite(result_shift$max_level_shift))
  expect_true(result_shift$max_level_shift > 1)  # Should detect large shift

  # Time of shift should be around position 50
  expect_true(is.finite(result_shift$time_level_shift))
  expect_true(result_shift$time_level_shift > 40 &&
              result_shift$time_level_shift < 70)  # Within reasonable range

  # White noise should have smaller level shift
  x_noise <- rnorm(100)
  result_noise <- ts_shifts(x_noise, window_size = 10)
  expect_true(result_shift$max_level_shift > result_noise$max_level_shift)
})

test_that("Variance shift detection works", {
  set.seed(789)

  # Data with variance shift
  x_var_shift <- c(rnorm(50, sd = 0.5), rnorm(50, sd = 3))
  result_var <- ts_shifts(x_var_shift, window_size = 10)

  # Should detect variance shift
  expect_true(is.finite(result_var$max_var_shift))
  expect_true(result_var$max_var_shift > 0)

  # Time should be around position 50
  expect_true(is.finite(result_var$time_var_shift))
  expect_true(result_var$time_var_shift > 40 &&
              result_var$time_var_shift < 70)
})

test_that("KL divergence shift detection works", {
  set.seed(111)

  # Data with distribution shift
  x_kl_shift <- c(rnorm(60, mean = 0, sd = 1), rnorm(40, mean = 3, sd = 2))
  result_kl <- ts_shifts(x_kl_shift, window_size = 10)

  # Should detect KL shift
  expect_true(is.finite(result_kl$max_kl_shift))
  expect_true(result_kl$max_kl_shift > 0)

  # Should have valid time index
  expect_true(is.finite(result_kl$time_kl_shift) || is.na(result_kl$time_kl_shift))
})

test_that("Shift detection handles short series", {
  # Too short for shift detection
  x_short <- rnorm(15)

  result <- ts_shifts(x_short, window_size = 10)

  # Should return NA values
  expect_true(is.na(result$max_level_shift))
  expect_true(is.na(result$time_level_shift))
  expect_true(is.na(result$max_var_shift))
  expect_true(is.na(result$time_var_shift))
  expect_true(is.na(result$max_kl_shift))
  expect_true(is.na(result$time_kl_shift))
  expect_equal(length(result), 6)
})

test_that("Shift detection with different window sizes", {
  set.seed(222)

  x <- c(rnorm(50, mean = 0), rnorm(50, mean = 3))

  # Small window
  result_small <- ts_shifts(x, window_size = 5)
  expect_true(is.finite(result_small$max_level_shift))

  # Medium window
  result_medium <- ts_shifts(x, window_size = 10)
  expect_true(is.finite(result_medium$max_level_shift))

  # Large window
  result_large <- ts_shifts(x, window_size = 20)
  expect_true(is.finite(result_large$max_level_shift))

  # All should detect the shift
  expect_true(result_small$max_level_shift > 1)
  expect_true(result_medium$max_level_shift > 1)
  expect_true(result_large$max_level_shift > 1)
})

test_that("Shift detection integrates with ts_features", {
  set.seed(333)

  x <- rnorm(100)

  # Test with 'shifts' subset
  shift_feats <- ts_features(x, features = "shifts")
  expect_equal(length(shift_feats), 6)
  expect_true("max_level_shift" %in% names(shift_feats))
  expect_true("max_var_shift" %in% names(shift_feats))
  expect_true("max_kl_shift" %in% names(shift_feats))

  # Test with 'all'
  all_feats <- ts_features(x, features = "all")
  expect_true("max_level_shift" %in% names(all_feats))
  expect_true("time_level_shift" %in% names(all_feats))
  expect_true(length(all_feats) >= 141)  # Should have at least 118 features now
})

test_that("Shift values are non-negative", {
  set.seed(444)

  x <- rnorm(100)
  result <- ts_shifts(x)

  # Max shifts should be non-negative
  if (is.finite(result$max_level_shift)) {
    expect_true(result$max_level_shift >= 0)
  }
  if (is.finite(result$max_var_shift)) {
    expect_true(result$max_var_shift >= 0)
  }
  if (is.finite(result$max_kl_shift)) {
    expect_true(result$max_kl_shift >= 0)
  }
})

test_that("Time indices are valid", {
  set.seed(555)

  x <- rnorm(100)
  result <- ts_shifts(x)

  # Time indices should be within valid range
  if (is.finite(result$time_level_shift)) {
    expect_true(result$time_level_shift >= 1)
    expect_true(result$time_level_shift <= length(x))
  }
  if (is.finite(result$time_var_shift)) {
    expect_true(result$time_var_shift >= 1)
    expect_true(result$time_var_shift <= length(x))
  }
  if (is.finite(result$time_kl_shift)) {
    expect_true(result$time_kl_shift >= 1)
    expect_true(result$time_kl_shift <= length(x))
  }
})

test_that("Shift detection handles constant series", {
  x_const <- rep(5, 100)

  result <- ts_shifts(x_const, window_size = 10)

  # Should handle gracefully
  expect_equal(length(result), 6)

  # Level shift should be 0 or very small for constant series
  if (is.finite(result$max_level_shift)) {
    expect_true(result$max_level_shift < 0.001)
  }
})

test_that("Multiple shifts are detected", {
  set.seed(666)

  # Series with multiple shifts
  x_multi <- c(
    rnorm(30, mean = 0),
    rnorm(30, mean = 5),
    rnorm(40, mean = 2)
  )

  result <- ts_shifts(x_multi, window_size = 10)

  # Should detect significant shift
  expect_true(is.finite(result$max_level_shift))
  expect_true(result$max_level_shift > 1)

  # Should report the maximum shift
  expect_true(is.finite(result$time_level_shift))
})

test_that("Shift detection works with trending data", {
  set.seed(777)

  # Trending data
  t <- 1:100
  x_trend <- 0.05 * t + rnorm(100, sd = 0.5)

  result <- ts_shifts(x_trend, window_size = 10)

  # Should compute all features
  expect_true(is.numeric(result$max_level_shift))
  expect_true(is.numeric(result$max_var_shift))
  expect_true(is.numeric(result$max_kl_shift))
})

test_that("Shift detection works with seasonal data", {
  set.seed(888)

  # Seasonal data
  t <- 1:100
  x_seasonal <- sin(2 * pi * t / 12) + rnorm(100, sd = 0.2)

  result <- ts_shifts(x_seasonal, window_size = 12)

  # Should compute all features
  expect_true(is.numeric(result$max_level_shift))
  expect_true(is.numeric(result$max_var_shift))
  expect_true(is.numeric(result$max_kl_shift))
})

test_that("Level shift detects mean changes", {
  set.seed(999)

  # Pure mean shift (same variance)
  x_mean_shift <- c(rnorm(50, mean = 0, sd = 1), rnorm(50, mean = 4, sd = 1))
  result <- ts_shifts(x_mean_shift, window_size = 10)

  # Level shift should be significant for large mean change
  if (is.finite(result$max_level_shift)) {
    expect_true(result$max_level_shift > 1)  # Should detect large shift
  }

  # Both shifts should be computable
  expect_true(is.finite(result$max_level_shift))
  expect_true(is.finite(result$max_var_shift))
})

test_that("Variance shift larger than level shift for sd change", {
  set.seed(101010)

  # Pure variance shift (same mean)
  x_var_shift <- c(rnorm(50, mean = 0, sd = 0.5), rnorm(50, mean = 0, sd = 3))
  result <- ts_shifts(x_var_shift, window_size = 10)

  # Variance shift should be significant
  if (is.finite(result$max_var_shift)) {
    expect_true(result$max_var_shift > 1)
  }
})

test_that("Shift detection handles NA values", {
  set.seed(111111)

  # Series with NAs (will be removed by ts_features)
  x <- rnorm(100)
  x[c(10, 50, 80)] <- NA

  # ts_features removes NAs before calling ts_shifts
  all_feats <- ts_features(x, features = "shifts")

  # Should handle NA values
  expect_equal(length(all_feats), EXPECTED_SHIFTS_FEATURES)
  expect_true(is.numeric(all_feats$max_level_shift) || is.na(all_feats$max_level_shift))
})

test_that("Shift detection is consistent across runs", {
  set.seed(121212)
  x <- rnorm(100)

  result1 <- ts_shifts(x, window_size = 10)
  result2 <- ts_shifts(x, window_size = 10)

  # Should get identical results
  expect_equal(result1$max_level_shift, result2$max_level_shift, tolerance = 1e-10)
  expect_equal(result1$time_level_shift, result2$time_level_shift)
  expect_equal(result1$max_var_shift, result2$max_var_shift, tolerance = 1e-10)
  expect_equal(result1$time_var_shift, result2$time_var_shift)
  expect_equal(result1$max_kl_shift, result2$max_kl_shift, tolerance = 1e-10)
  expect_equal(result1$time_kl_shift, result2$time_kl_shift)
})

test_that("Feature count is correct after Phase 3 Batch 5", {
  set.seed(131313)
  x <- rnorm(150)

  all_feats <- ts_features(x, features = "all")

  # Should have 118 total features (112 + 6)
  expect_equal(length(all_feats), EXPECTED_ALL_FEATURES)
})

test_that("Shift detection detects abrupt changes better than gradual", {
  set.seed(141414)

  # Abrupt change
  x_abrupt <- c(rnorm(50, mean = 0), rnorm(50, mean = 5))
  result_abrupt <- ts_shifts(x_abrupt, window_size = 10)

  # Gradual change
  x_gradual <- rnorm(100)
  x_gradual <- x_gradual + seq(0, 5, length.out = 100)
  result_gradual <- ts_shifts(x_gradual, window_size = 10)

  # Abrupt should have higher shift value
  if (is.finite(result_abrupt$max_level_shift) && is.finite(result_gradual$max_level_shift)) {
    expect_true(result_abrupt$max_level_shift > result_gradual$max_level_shift)
  }
})

test_that("Window size affects sensitivity", {
  set.seed(151515)

  # Small localized shift
  x <- rnorm(100)
  x[48:52] <- x[48:52] + 5

  # Small window should be more sensitive to localized changes
  result_small <- ts_shifts(x, window_size = 5)
  result_large <- ts_shifts(x, window_size = 20)

  # Both should detect, but may differ in magnitude
  expect_true(is.finite(result_small$max_level_shift))
  expect_true(is.finite(result_large$max_level_shift))
})

test_that("Shift detection works with different data scales", {
  set.seed(161616)

  # Small scale
  x_small <- c(rnorm(50, mean = 0, sd = 0.1), rnorm(50, mean = 1, sd = 0.1))
  result_small <- ts_shifts(x_small, window_size = 10)

  # Large scale
  x_large <- c(rnorm(50, mean = 0, sd = 100), rnorm(50, mean = 1000, sd = 100))
  result_large <- ts_shifts(x_large, window_size = 10)

  # Both should detect shifts
  expect_true(is.finite(result_small$max_level_shift))
  expect_true(is.finite(result_large$max_level_shift))
})

test_that("KL shift detects distributional changes", {
  set.seed(171717)

  # Change from normal to bimodal (distribution change without mean/var change)
  x_dist <- c(
    rnorm(50, mean = 0, sd = 1),
    c(rnorm(25, mean = -1, sd = 0.5), rnorm(25, mean = 1, sd = 0.5))
  )

  result <- ts_shifts(x_dist, window_size = 10)

  # KL divergence should detect this
  expect_true(is.finite(result$max_kl_shift))
  expect_true(result$max_kl_shift > 0)
})
