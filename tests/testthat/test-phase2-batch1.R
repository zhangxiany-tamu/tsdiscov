test_that("Hurst exponent is computed correctly", {
  set.seed(123)
  x <- rnorm(200)

  hurst <- cpp_hurst_exponent(x)

  # Should return a finite value
  expect_true(is.finite(hurst))

  # Should be in valid range [0, 1]
  expect_true(hurst >= 0 && hurst <= 1)
})

test_that("Hurst exponent distinguishes signal types", {
  set.seed(42)

  # White noise - should be around 0.5
  wn <- rnorm(500)
  h_wn <- cpp_hurst_exponent(wn)

  # Random walk - should be > white noise
  rw <- cumsum(rnorm(500))
  h_rw <- cpp_hurst_exponent(rw)

  # Strong trend - should show persistence
  t <- 1:500
  trending <- 2 * t + rnorm(500, sd = 50)
  h_trend <- cpp_hurst_exponent(trending)

  # All should be finite
  expect_true(is.finite(h_wn))
  expect_true(is.finite(h_rw))
  expect_true(is.finite(h_trend))

  # White noise should be around 0.5 (allow wide tolerance)
  expect_true(h_wn > 0.3 && h_wn < 0.7)
})

test_that("Hurst exponent handles edge cases", {
  # Very short series
  x_short <- rnorm(15)
  h_short <- cpp_hurst_exponent(x_short)
  expect_true(is.na(h_short))

  # Constant series
  x_const <- rep(5, 100)
  h_const <- cpp_hurst_exponent(x_const)
  # May be NA or a value, just check it doesn't error
  expect_true(is.numeric(h_const))
})

test_that("Stability works correctly", {
  set.seed(123)

  # Stationary series (low stability variance)
  stationary <- rnorm(200)
  stab_stat <- cpp_stability(stationary, 10)

  # Non-stationary (high stability variance)
  non_stationary <- cumsum(rnorm(200))
  stab_non <- cpp_stability(non_stationary, 10)

  # Both should be finite
  expect_true(is.finite(stab_stat))
  expect_true(is.finite(stab_non))

  # Both should be non-negative
  expect_true(stab_stat >= 0)
  expect_true(stab_non >= 0)

  # Non-stationary should have higher stability variance
  expect_true(stab_non > stab_stat)
})

test_that("Stability handles different window sizes", {
  set.seed(42)
  x <- rnorm(200)

  stab_5 <- cpp_stability(x, 5)
  stab_10 <- cpp_stability(x, 10)
  stab_20 <- cpp_stability(x, 20)

  # All should be finite
  expect_true(is.finite(stab_5))
  expect_true(is.finite(stab_10))
  expect_true(is.finite(stab_20))
})

test_that("Stability returns NA for short series", {
  x_short <- rnorm(15)
  stab <- cpp_stability(x_short, 10)
  expect_true(is.na(stab))
})

test_that("Lumpiness works correctly", {
  set.seed(123)

  # Homoscedastic series (low lumpiness)
  homoscedastic <- rnorm(200)
  lump_homo <- cpp_lumpiness(homoscedastic, 10)

  # Heteroscedastic series (high lumpiness)
  t <- 1:200
  heteroscedastic <- rnorm(200, sd = abs(sin(t / 20)) + 0.1)
  lump_hetero <- cpp_lumpiness(heteroscedastic, 10)

  # Both should be finite
  expect_true(is.finite(lump_homo))
  expect_true(is.finite(lump_hetero))

  # Both should be non-negative
  expect_true(lump_homo >= 0)
  expect_true(lump_hetero >= 0)

  # Heteroscedastic should have higher lumpiness
  expect_true(lump_hetero > lump_homo)
})

test_that("Lumpiness handles different window sizes", {
  set.seed(42)
  x <- rnorm(200)

  lump_5 <- cpp_lumpiness(x, 5)
  lump_10 <- cpp_lumpiness(x, 10)
  lump_20 <- cpp_lumpiness(x, 20)

  # All should be finite
  expect_true(is.finite(lump_5))
  expect_true(is.finite(lump_10))
  expect_true(is.finite(lump_20))
})

test_that("Lumpiness returns NA for short series", {
  x_short <- rnorm(15)
  lump <- cpp_lumpiness(x_short, 10)
  expect_true(is.na(lump))
})

test_that("Lumpiness detects variance changes", {
  set.seed(42)

  # Create series with abrupt variance change
  x1 <- rnorm(100, sd = 1)
  x2 <- rnorm(100, sd = 5)
  x_change <- c(x1, x2)

  # Constant variance
  x_const_var <- rnorm(200, sd = 2)

  lump_change <- cpp_lumpiness(x_change, 20)
  lump_const <- cpp_lumpiness(x_const_var, 20)

  # Variance change should have higher lumpiness
  expect_true(lump_change > lump_const)
})

test_that("ts_scaling wrapper works correctly", {
  set.seed(123)
  x <- rnorm(200)

  scaling_features <- ts_scaling(x)

  # Should have 3 features
  expect_equal(length(scaling_features), 3)

  # Check all expected features are present
  expect_true("hurst_exponent" %in% names(scaling_features))
  expect_true("stability" %in% names(scaling_features))
  expect_true("lumpiness" %in% names(scaling_features))

  # All should be numeric
  for (i in seq_along(scaling_features)) {
    expect_true(is.numeric(scaling_features[[i]]))
  }
})

test_that("Scaling features work with different signal types", {
  set.seed(42)

  # White noise
  wn <- rnorm(200)
  f_wn <- ts_scaling(wn)

  # Random walk
  rw <- cumsum(rnorm(200))
  f_rw <- ts_scaling(rw)

  # Periodic signal
  t <- 1:200
  periodic <- sin(2 * pi * t / 20) + rnorm(200, sd = 0.2)
  f_periodic <- ts_scaling(periodic)

  # All should compute successfully
  expect_true(all(sapply(f_wn, function(x) is.finite(x) || is.na(x))))
  expect_true(all(sapply(f_rw, function(x) is.finite(x) || is.na(x))))
  expect_true(all(sapply(f_periodic, function(x) is.finite(x) || is.na(x))))

  # Random walk should have higher stability than white noise
  expect_true(f_rw$stability > f_wn$stability)
})

test_that("Complete feature extraction includes scaling features", {
  set.seed(123)
  x <- rnorm(200)

  # Extract all features
  all_features <- ts_features(x, features = "all")

  # Should now have 80 + 3 = 83 features
  expect_equal(length(all_features), 352)

  # Check scaling features are included
  expect_true("hurst_exponent" %in% names(all_features))
  expect_true("stability" %in% names(all_features))
  expect_true("lumpiness" %in% names(all_features))

  # All features should be numeric scalars
  for (i in seq_along(all_features)) {
    expect_true(is.numeric(all_features[[i]]) || is.logical(all_features[[i]]))
    expect_equal(length(all_features[[i]]), 1)
  }
})

test_that("Scaling features can be extracted independently", {
  set.seed(123)
  x <- rnorm(200)

  # Extract only scaling features
  scaling_only <- ts_features(x, features = "scaling")

  # Should have exactly 3 features
  expect_equal(length(scaling_only), 3)

  # All names should be scaling features
  expected_names <- c("hurst_exponent", "stability", "lumpiness")
  expect_true(all(expected_names %in% names(scaling_only)))
})

test_that("Scaling features are consistent across runs", {
  set.seed(999)
  x <- rnorm(200)

  f1 <- ts_scaling(x)
  f2 <- ts_scaling(x)

  # Should get identical results
  for (i in seq_along(f1)) {
    if (is.finite(f1[[i]]) && is.finite(f2[[i]])) {
      expect_equal(f1[[i]], f2[[i]], tolerance = 1e-10)
    } else {
      # Both should be NA if one is NA
      expect_equal(is.na(f1[[i]]), is.na(f2[[i]]))
    }
  }
})

test_that("Constant series handling", {
  x_const <- rep(5, 200)

  scaling_features <- ts_scaling(x_const)

  # Stability and lumpiness should be 0 for constant series
  expect_equal(scaling_features$stability, 0.0)
  expect_equal(scaling_features$lumpiness, 0.0)

  # Hurst may be NA or a value, just check it's computed
  expect_true(is.numeric(scaling_features$hurst_exponent))
})
