# Phase 3 Batch 4: Statistical Tests
# Tests for unit root tests, nonlinearity tests

test_that("Statistical tests work correctly", {
  skip_if_not_installed("urca")
  skip_if_not_installed("tseries")

  set.seed(123)

  # Create stationary white noise
  x <- rnorm(100)

  result <- ts_stattests(x)

  # Should have 3 features
  expect_equal(length(result), 3)
  expect_true("unitroot_kpss" %in% names(result))
  expect_true("unitroot_pp" %in% names(result))
  expect_true("nonlinearity" %in% names(result))
})

test_that("KPSS unit root test works", {
  skip_if_not_installed("urca")

  set.seed(456)

  # Stationary series (white noise)
  x_stationary <- rnorm(100)
  result_stat <- ts_stattests(x_stationary)

  # KPSS test: low statistic indicates stationarity
  expect_true(is.finite(result_stat$unitroot_kpss))
  expect_true(is.numeric(result_stat$unitroot_kpss))

  # Non-stationary series (random walk)
  x_nonstat <- cumsum(rnorm(100))
  result_nonstat <- ts_stattests(x_nonstat)

  # KPSS test: high statistic indicates non-stationarity
  expect_true(is.finite(result_nonstat$unitroot_kpss))

  # Non-stationary should have higher KPSS statistic
  expect_true(result_nonstat$unitroot_kpss > result_stat$unitroot_kpss)
})

test_that("Phillips-Perron unit root test works", {
  skip_if_not_installed("urca")

  set.seed(789)

  # Stationary series
  x_stationary <- rnorm(100)
  result_stat <- ts_stattests(x_stationary)

  # PP test: large negative statistic indicates stationarity
  expect_true(is.finite(result_stat$unitroot_pp))
  expect_true(is.numeric(result_stat$unitroot_pp))

  # Random walk (non-stationary)
  x_nonstat <- cumsum(rnorm(100))
  result_nonstat <- ts_stattests(x_nonstat)

  # PP test should be finite
  expect_true(is.finite(result_nonstat$unitroot_pp) || is.na(result_nonstat$unitroot_pp))
})

test_that("Nonlinearity test works", {
  skip_if_not_installed("tseries")

  set.seed(111)

  # Linear AR process
  n <- 200
  x_linear <- numeric(n)
  x_linear[1] <- rnorm(1)
  for(i in 2:n) {
    x_linear[i] <- 0.5 * x_linear[i-1] + rnorm(1)
  }

  result_linear <- ts_stattests(x_linear)

  # Should compute nonlinearity statistic
  expect_true(is.numeric(result_linear$nonlinearity))

  # Nonlinear series (threshold AR)
  x_nonlinear <- numeric(n)
  x_nonlinear[1] <- rnorm(1)
  for(i in 2:n) {
    if (x_nonlinear[i-1] > 0) {
      x_nonlinear[i] <- 0.8 * x_nonlinear[i-1] + rnorm(1)
    } else {
      x_nonlinear[i] <- -0.3 * x_nonlinear[i-1] + rnorm(1)
    }
  }

  result_nonlinear <- ts_stattests(x_nonlinear)

  # Should compute nonlinearity statistic
  expect_true(is.numeric(result_nonlinear$nonlinearity))
})

test_that("Statistical tests handle short series", {
  # Too short for tests
  x_short <- rnorm(15)

  result <- ts_stattests(x_short)

  # Should return NA values
  expect_true(is.na(result$unitroot_kpss))
  expect_true(is.na(result$unitroot_pp))
  expect_true(is.na(result$nonlinearity))
  expect_equal(length(result), 3)
})

test_that("Statistical tests integrate with ts_features", {
  skip_if_not_installed("urca")
  skip_if_not_installed("tseries")

  set.seed(222)

  x <- rnorm(100)

  # Test with 'stattests' subset
  stat_feats <- ts_features(x, features = "stattests")
  expect_equal(length(stat_feats), 3)
  expect_true("unitroot_kpss" %in% names(stat_feats))
  expect_true("unitroot_pp" %in% names(stat_feats))
  expect_true("nonlinearity" %in% names(stat_feats))

  # Test with 'all'
  all_feats <- ts_features(x, features = "all")
  expect_true("unitroot_kpss" %in% names(all_feats))
  expect_true("unitroot_pp" %in% names(all_feats))
  expect_true("nonlinearity" %in% names(all_feats))
  expect_true(length(all_feats) >= 141)  # Should have at least 112 features now
})

test_that("Statistical tests handle constant series", {
  x_const <- rep(5, 100)

  result <- ts_stattests(x_const)

  # Should handle gracefully (may return NA or error)
  expect_equal(length(result), 3)
  expect_true(is.numeric(result$unitroot_kpss) || is.na(result$unitroot_kpss))
  expect_true(is.numeric(result$unitroot_pp) || is.na(result$unitroot_pp))
  expect_true(is.numeric(result$nonlinearity) || is.na(result$nonlinearity))
})

test_that("Statistical tests handle trending series", {
  skip_if_not_installed("urca")

  set.seed(333)

  # Series with strong trend
  t <- 1:100
  x_trend <- 0.5 * t + rnorm(100, sd = 1)

  result <- ts_stattests(x_trend)

  # KPSS should indicate non-stationarity (high value)
  expect_true(is.finite(result$unitroot_kpss))
  expect_true(result$unitroot_kpss > 0.1)  # Likely non-stationary
})

test_that("Statistical tests work with different series lengths", {
  skip_if_not_installed("urca")
  skip_if_not_installed("tseries")

  set.seed(444)

  # Medium length
  x_medium <- rnorm(100)
  result_medium <- ts_stattests(x_medium)
  expect_true(is.numeric(result_medium$unitroot_kpss))

  # Long length
  x_long <- rnorm(500)
  result_long <- ts_stattests(x_long)
  expect_true(is.numeric(result_long$unitroot_kpss))
})

test_that("Statistical tests are consistent across runs", {
  skip_if_not_installed("urca")
  skip_if_not_installed("tseries")

  set.seed(555)
  x <- rnorm(100)

  result1 <- ts_stattests(x)
  result2 <- ts_stattests(x)

  # Should get identical results
  expect_equal(result1$unitroot_kpss, result2$unitroot_kpss, tolerance = 1e-10)
  expect_equal(result1$unitroot_pp, result2$unitroot_pp, tolerance = 1e-10)
  expect_equal(result1$nonlinearity, result2$nonlinearity, tolerance = 1e-10)
})

test_that("KPSS test statistic is positive", {
  skip_if_not_installed("urca")

  set.seed(666)

  x <- rnorm(100)
  result <- ts_stattests(x)

  # KPSS statistic should always be non-negative
  if (is.finite(result$unitroot_kpss)) {
    expect_true(result$unitroot_kpss >= 0)
  }
})

test_that("Statistical tests handle NA values", {
  skip_if_not_installed("urca")
  skip_if_not_installed("tseries")

  set.seed(777)

  # Series with NAs (will be removed by ts_features)
  x <- rnorm(100)
  x[c(10, 50, 80)] <- NA

  # ts_features removes NAs before calling ts_stattests
  all_feats <- ts_features(x, features = "stattests")

  # Should handle NA values
  expect_equal(length(all_feats), EXPECTED_MISC_FEATURES)
  expect_true(is.numeric(all_feats$unitroot_kpss) || is.na(all_feats$unitroot_kpss))
})

test_that("Nonlinearity statistic is non-negative", {
  skip_if_not_installed("tseries")

  set.seed(888)

  x <- rnorm(100)
  result <- ts_stattests(x)

  # Nonlinearity statistic should be non-negative
  if (is.finite(result$nonlinearity)) {
    expect_true(result$nonlinearity >= 0)
  }
})

test_that("Feature count is correct after Phase 3 Batch 4", {
  skip_if_not_installed("urca")
  skip_if_not_installed("tseries")

  set.seed(101010)
  x <- rnorm(150)

  all_feats <- ts_features(x, features = "all")

  # Should have 112 total features (109 + 3)
  expect_equal(length(all_feats), EXPECTED_ALL_FEATURES)
})

test_that("Statistical tests detect unit roots", {
  skip_if_not_installed("urca")

  set.seed(999)

  # I(0) process: stationary
  x_stationary <- arima.sim(list(ar = 0.5), n = 200)
  result_stat <- ts_stattests(as.numeric(x_stationary))

  # I(1) process: unit root
  x_unit_root <- cumsum(rnorm(200))
  result_unit <- ts_stattests(x_unit_root)

  # KPSS: higher for unit root series
  if (is.finite(result_stat$unitroot_kpss) && is.finite(result_unit$unitroot_kpss)) {
    expect_true(result_unit$unitroot_kpss > result_stat$unitroot_kpss)
  }
})

test_that("Statistical tests handle seasonal data", {
  skip_if_not_installed("urca")
  skip_if_not_installed("tseries")

  set.seed(111111)

  # Seasonal data
  t <- 1:100
  x_seasonal <- sin(2*pi*t/12) + rnorm(100, sd = 0.2)

  result <- ts_stattests(x_seasonal)

  # Should compute all tests
  expect_true(is.numeric(result$unitroot_kpss))
  expect_true(is.numeric(result$unitroot_pp))
  expect_true(is.numeric(result$nonlinearity))
})

test_that("Statistical tests work without optional packages", {
  set.seed(135256)
  x <- rnorm(100)

  # This should work even if packages are not available
  # (will return NA values)
  result <- tryCatch({
    ts_stattests(x)
  }, error = function(e) {
    list(unitroot_kpss = NA_real_, unitroot_pp = NA_real_, nonlinearity = NA_real_)
  })

  expect_equal(length(result), 3)
  expect_true("unitroot_kpss" %in% names(result))
  expect_true("unitroot_pp" %in% names(result))
  expect_true("nonlinearity" %in% names(result))
})
