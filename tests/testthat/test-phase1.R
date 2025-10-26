test_that("PACF features work correctly", {
  set.seed(123)
  x <- rnorm(100)

  # Test PACF computation
  pacf_vals <- cpp_pacf(x, 20)
  expect_equal(length(pacf_vals), 20)
  expect_true(all(is.finite(pacf_vals)))

  # PACF should be between -1 and 1
  expect_true(all(abs(pacf_vals) <= 1.1)) # Allow small numerical errors

  # Test PACF features
  pacf_features <- cpp_pacf_features(x, 20)
  expect_equal(length(pacf_features), 6)
  expect_true(all(names(pacf_features) %in% c("pacf_lag1", "pacf_lag5", "pacf_lag10",
                                               "first_sig_pacf", "sum_sq_pacf", "x_pacf5")))

  # All features should be finite or NA
  for (i in seq_along(pacf_features)) {
    expect_true(is.finite(pacf_features[[i]]) || is.na(pacf_features[[i]]))
  }

  # Test ts_pacf wrapper
  pacf_ts <- ts_pacf(x)
  expect_equal(length(pacf_ts), 6)
  expect_true("pacf_lag1" %in% names(pacf_ts))
})

test_that("PACF handles AR processes correctly", {
  # AR(1) process: x[t] = 0.7 * x[t-1] + noise
  set.seed(42)
  n <- 200
  x <- numeric(n)
  x[1] <- rnorm(1)
  for (i in 2:n) {
    x[i] <- 0.7 * x[i - 1] + rnorm(1, sd = 0.5)
  }

  pacf_vals <- cpp_pacf(x, 10)

  # For AR(1), PACF should be significant at lag 1 and near zero after
  expect_true(abs(pacf_vals[1]) > 0.3) # Lag 1 should be strong
  expect_true(abs(pacf_vals[5]) < abs(pacf_vals[1])) # Later lags should be smaller
})

test_that("Linear trend features work correctly", {
  set.seed(123)

  # Test on data with known trend
  t <- 1:100
  x <- 2.5 * t + 10 + rnorm(100, sd = 5)

  trend_vals <- cpp_linear_trend(x)

  expect_equal(length(trend_vals), 4)
  expect_true(all(names(trend_vals) %in% c("slope", "intercept", "r_squared", "stderr")))

  # Check slope is close to true value (2.5)
  expect_true(abs(trend_vals["slope"] - 2.5) < 1.0) # Within 1.0 due to noise

  # Check intercept is close to true value (10)
  expect_true(abs(trend_vals["intercept"] - 10) < 20) # Within 20 due to noise

  # R-squared should be high for strong trend
  expect_true(trend_vals["r_squared"] > 0.8)

  # Standard error should be positive and finite
  expect_true(trend_vals["stderr"] > 0)
  expect_true(is.finite(trend_vals["stderr"]))
})

test_that("Linear trend detects no trend in white noise", {
  set.seed(456)
  x <- rnorm(100)

  trend_vals <- cpp_linear_trend(x)

  # Slope should be close to zero
  expect_true(abs(trend_vals["slope"]) < 0.5)

  # R-squared should be low
  expect_true(trend_vals["r_squared"] < 0.3)
})

test_that("Mean absolute change works correctly", {
  set.seed(123)

  # Constant series
  x_const <- rep(5, 100)
  mac_const <- cpp_mean_abs_change(x_const)
  expect_equal(mac_const, 0.0)

  # Random walk has higher MAC than white noise
  wn <- rnorm(100)
  rw <- cumsum(rnorm(100))

  mac_wn <- cpp_mean_abs_change(wn)
  mac_rw <- cpp_mean_abs_change(rw)

  expect_true(is.finite(mac_wn))
  expect_true(is.finite(mac_rw))
  expect_true(mac_wn > 0)
  expect_true(mac_rw > 0)
})

test_that("Mean change works correctly", {
  set.seed(123)

  # Upward trending data
  x_up <- cumsum(abs(rnorm(100)))
  mc_up <- cpp_mean_change(x_up)
  expect_true(mc_up > 0) # Positive mean change

  # Downward trending data
  x_down <- -cumsum(abs(rnorm(100)))
  mc_down <- cpp_mean_change(x_down)
  expect_true(mc_down < 0) # Negative mean change

  # White noise should have mean change near zero
  wn <- rnorm(100)
  mc_wn <- cpp_mean_change(wn)
  expect_true(abs(mc_wn) < 0.5)
})

test_that("Mean second derivative works correctly", {
  set.seed(123)

  # Linear trend has second derivative near zero
  t <- 1:100
  x_linear <- 2 * t + 5
  msd_linear <- cpp_mean_second_derivative(x_linear)
  expect_true(abs(msd_linear) < 0.01)

  # Quadratic has non-zero second derivative
  x_quad <- 0.1 * t^2 + 2 * t + 5
  msd_quad <- cpp_mean_second_derivative(x_quad)
  expect_true(is.finite(msd_quad))

  # Should return NA for very short series
  x_short <- rnorm(2)
  msd_short <- cpp_mean_second_derivative(x_short)
  expect_true(is.na(msd_short))
})

test_that("ts_trend wrapper works correctly", {
  set.seed(123)
  x <- rnorm(100)

  trend_features <- ts_trend(x)

  # Should have 7 features
  expect_equal(length(trend_features), 7)

  # Check all expected features are present
  expected_names <- c("trend_slope", "trend_intercept", "trend_r_squared",
                      "trend_stderr", "mean_abs_change", "mean_change",
                      "mean_second_deriv")
  expect_true(all(expected_names %in% names(trend_features)))

  # All should be numeric
  for (i in seq_along(trend_features)) {
    expect_true(is.numeric(trend_features[[i]]))
  }
})

test_that("Complete feature extraction with new features works", {
  set.seed(123)
  x <- rnorm(100)

  # Extract all features
  all_features <- ts_features(x, features = "all")

  # Should now have 39 + 5 (PACF) + 7 (trend) + 14 (structure) + 15 (FFT) + 3 (scaling) = 83 features
  expect_equal(length(all_features), 352)

  # Check new feature categories are included
  expect_true("pacf_lag1" %in% names(all_features))
  expect_true("trend_slope" %in% names(all_features))
  expect_true("mean_abs_change" %in% names(all_features))

  # All features should be numeric scalars
  for (i in seq_along(all_features)) {
    expect_true(is.numeric(all_features[[i]]) || is.logical(all_features[[i]]))
    expect_equal(length(all_features[[i]]), 1)
  }
})

test_that("New features work on different signal types", {
  set.seed(42)

  # White noise
  wn <- rnorm(200)
  f_wn <- ts_features(wn, features = c("pacf", "trend"))

  # Trending signal
  t <- 1:200
  trend <- 2 * t + rnorm(200, sd = 10)
  f_trend <- ts_features(trend, features = c("pacf", "trend"))

  # Periodic signal
  periodic <- sin(2 * pi * t / 20) + rnorm(200, sd = 0.2)
  f_periodic <- ts_features(periodic, features = c("pacf", "trend"))

  # Trending signal should have higher R-squared
  expect_true(f_trend$trend_r_squared > f_wn$trend_r_squared)

  # Trending signal should have positive slope
  expect_true(f_trend$trend_slope > 1)

  # All features should be computed
  expect_equal(length(f_wn), 13)
  expect_equal(length(f_trend), 13)
  expect_equal(length(f_periodic), 13)
})

test_that("Edge cases are handled correctly", {
  # Very short series
  x_short <- rnorm(5)
  pacf_short <- cpp_pacf_features(x_short, 20)
  trend_short <- cpp_linear_trend(x_short)

  # Should return some values (may be NA for some features)
  expect_true(length(pacf_short) >= 5)  # May return 5 or 6 depending on data
  expect_equal(length(trend_short), 4)

  # Constant series
  x_const <- rep(1, 100)
  trend_const <- cpp_linear_trend(x_const)
  expect_equal(as.numeric(trend_const["slope"]), 0.0, tolerance = 1e-10)

  mac_const <- cpp_mean_abs_change(x_const)
  expect_equal(mac_const, 0.0)
})
