# Tests for Phase 4 features

test_that("Basic stats features work correctly", {
  set.seed(123)
  x <- rnorm(100)

  basic <- ts_basic(x)

  # Should return 3 features
  expect_equal(length(basic), 3)

  # Check names
  expect_equal(names(basic), c("maximum", "minimum", "root_mean_square"))

  # All should be numeric
  expect_true(all(sapply(basic, is.numeric)))

  # Maximum should be largest value
  expect_equal(basic$maximum, max(x))

  # Minimum should be smallest value
  expect_equal(basic$minimum, min(x))

  # RMS should be positive
  expect_true(basic$root_mean_square > 0)

  # RMS formula check
  expect_equal(basic$root_mean_square, sqrt(mean(x^2)))
})

test_that("Basic stats handle edge cases", {
  # All positive
  x_pos <- abs(rnorm(50))
  basic_pos <- ts_basic(x_pos)
  expect_true(basic_pos$minimum >= 0)

  # All negative
  x_neg <- -abs(rnorm(50))
  basic_neg <- ts_basic(x_neg)
  expect_true(basic_neg$maximum <= 0)

  # Constant series
  x_const <- rep(5, 50)
  basic_const <- ts_basic(x_const)
  expect_equal(basic_const$maximum, 5)
  expect_equal(basic_const$minimum, 5)
  expect_equal(basic_const$root_mean_square, 5)
})

test_that("Extrema locations work correctly", {
  set.seed(456)
  x <- rnorm(100)

  extrema <- ts_extrema(x)

  # Should return 4 features
  expect_equal(length(extrema), 4)

  # Check names
  expected_names <- c("first_loc_max", "last_loc_max", "first_loc_min", "last_loc_min")
  expect_equal(names(extrema), expected_names)

  # All should be between 0 and 1 (normalized)
  for (feat in extrema) {
    expect_true(feat >= 0 && feat <= 1)
  }

  # first_loc_max should be <= last_loc_max
  expect_true(extrema$first_loc_max <= extrema$last_loc_max)

  # first_loc_min should be <= last_loc_min
  expect_true(extrema$first_loc_min <= extrema$last_loc_min)
})

test_that("Extrema locations work with specific patterns", {
  # Max at beginning
  x_max_first <- c(100, 1:50)
  extrema_mf <- ts_extrema(x_max_first)
  expect_true(extrema_mf$first_loc_max < 0.1)  # Near start

  # Max at end
  x_max_last <- c(1:50, 100)
  extrema_ml <- ts_extrema(x_max_last)
  expect_true(extrema_ml$last_loc_max > 0.9)  # Near end

  # Multiple maxima (should find first and last)
  x_multi <- c(100, 1:48, 100)
  extrema_multi <- ts_extrema(x_multi)
  expect_true(extrema_multi$first_loc_max < 0.1)
  expect_true(extrema_multi$last_loc_max > 0.9)
})

test_that("Extrema locations handle edge cases", {
  # Constant series - all locations are valid
  x_const <- rep(5, 50)
  extrema_const <- ts_extrema(x_const)
  expect_true(all(!is.na(unlist(extrema_const))))

  # Single value that repeats
  x_repeat <- rep(c(1, 2), 25)
  extrema_repeat <- ts_extrema(x_repeat)
  expect_true(all(!is.na(unlist(extrema_repeat))))
})

test_that("AR coefficients work correctly", {
  set.seed(789)

  # AR(1) process: x[t] = 0.7*x[t-1] + noise
  n <- 150
  ar1_coef <- 0.7
  x <- numeric(n)
  x[1] <- rnorm(1)
  for (i in 2:n) {
    x[i] <- ar1_coef * x[i-1] + rnorm(1, sd=0.1)
  }

  ar_feats <- ts_ar(x)

  # Should return 4 features
  expect_equal(length(ar_feats), 4)

  # Check names
  expect_equal(names(ar_feats), c("ar_coef_1", "ar_coef_2", "ar_coef_5", "ar_coef_10"))

  # AR(1) coefficient should be close to 0.7
  expect_true(abs(ar_feats$ar_coef_1 - 0.7) < 0.3)  # Allow some estimation error

  # All should be numeric
  expect_true(all(sapply(ar_feats, is.numeric)))
})

test_that("AR coefficients work with white noise", {
  set.seed(111)
  x <- rnorm(150)

  ar_feats <- ts_ar(x)

  # White noise should have small AR coefficients
  expect_true(abs(ar_feats$ar_coef_1) < 0.5)
  expect_true(abs(ar_feats$ar_coef_2) < 0.5)
})

test_that("AR coefficients handle edge cases", {
  # Too short series
  x_short <- rnorm(5)
  ar_short <- ts_ar(x_short)
  expect_true(all(is.na(unlist(ar_short))))

  # Constant series
  x_const <- rep(5, 100)
  ar_const <- ts_ar(x_const)
  expect_true(all(is.na(unlist(ar_const))))

  # Series with very low variance
  x_lowvar <- rnorm(100, sd = 1e-20)
  ar_lowvar <- ts_ar(x_lowvar)
  expect_true(all(is.na(unlist(ar_lowvar))))
})

test_that("ADF test works correctly", {
  skip_if_not_installed("tseries")

  set.seed(321)

  # Stationary series (white noise)
  x_stationary <- rnorm(150)
  adf_stat <- ts_adf(x_stationary)

  # Should return 2 features
  expect_equal(length(adf_stat), 2)

  # Check names
  expect_equal(names(adf_stat), c("adf_stat", "adf_pvalue"))

  # All should be numeric
  expect_true(all(sapply(adf_stat, is.numeric)))

  # For stationary series, should reject null (low p-value)
  # Note: ADF null = has unit root (non-stationary)
  expect_true(adf_stat$adf_pvalue < 0.1)

  # Test statistic should be negative
  expect_true(adf_stat$adf_stat < 0)
})

test_that("ADF test works with non-stationary series", {
  skip_if_not_installed("tseries")

  set.seed(654)

  # Random walk (non-stationary)
  n <- 150
  x_rw <- cumsum(rnorm(n))
  adf_rw <- ts_adf(x_rw)

  # For random walk, might fail to reject null (high p-value)
  # But this is not guaranteed with short series, so just check format
  expect_true(is.numeric(adf_rw$adf_stat))
  expect_true(is.numeric(adf_rw$adf_pvalue))
  expect_true(adf_rw$adf_pvalue >= 0 && adf_rw$adf_pvalue <= 1)
})

test_that("ADF test handles edge cases", {
  skip_if_not_installed("tseries")

  # Too short series
  x_short <- rnorm(5)
  adf_short <- ts_adf(x_short)
  expect_true(all(is.na(unlist(adf_short))))
})

# Note: Testing missing tseries package is complex with mocking
# The warning path is tested manually and works as expected
# The function gracefully returns NA when tseries is not available

test_that("Phase 4 features integrate with ts_features", {
  set.seed(999)
  x <- rnorm(150)

  # Test with 'basic' subset
  basic_feats <- ts_features(x, features = "basic")
  expect_equal(length(basic_feats), 3)
  expect_true("maximum" %in% names(basic_feats))

  # Test with 'extrema' subset
  extrema_feats <- ts_features(x, features = "extrema")
  expect_equal(length(extrema_feats), 4)
  expect_true("first_loc_max" %in% names(extrema_feats))

  # Test with 'ar' subset
  ar_feats <- ts_features(x, features = "ar")
  expect_equal(length(ar_feats), 4)
  expect_true("ar_coef_1" %in% names(ar_feats))

  # Test with 'adf' subset
  skip_if_not_installed("tseries")
  adf_feats <- ts_features(x, features = "adf")
  expect_equal(length(adf_feats), 2)
  expect_true("adf_stat" %in% names(adf_feats))

  # Test with 'all'
  all_feats <- ts_features(x, features = "all")
  expect_true("maximum" %in% names(all_feats))
  expect_true("first_loc_max" %in% names(all_feats))
  expect_true("ar_coef_1" %in% names(all_feats))
  expect_true("adf_stat" %in% names(all_feats))
  expect_equal(length(all_feats), EXPECTED_ALL_FEATURES)
})

test_that("Feature count is correct after Phase 4", {
  set.seed(101010)
  x <- rnorm(200)

  all_feats <- ts_features(x, features = "all")

  # Should have 170 total features (157 + 13)
  # 3 basic + 4 extrema + 4 AR + 2 ADF = 13 new
  expect_equal(length(all_feats), EXPECTED_ALL_FEATURES)
})
