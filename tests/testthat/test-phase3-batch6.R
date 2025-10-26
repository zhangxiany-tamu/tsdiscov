# Phase 3 Batch 6: Miscellaneous Features
# Tests for firstzero_ac, zero_proportion, std1st_der

test_that("Miscellaneous features work correctly", {
  set.seed(123)

  x <- rnorm(100)

  result <- ts_misc(x)

  # Should have 3 features
  expect_equal(length(result), 3)
  expect_true("firstzero_ac" %in% names(result))
  expect_true("zero_proportion" %in% names(result))
  expect_true("std1st_der" %in% names(result))
})

test_that("First zero ACF crossing works", {
  set.seed(456)

  # White noise - ACF should cross zero quickly
  x_wn <- rnorm(100)
  result_wn <- ts_misc(x_wn)

  expect_true(is.numeric(result_wn$firstzero_ac))
  expect_true(result_wn$firstzero_ac >= 0)

  # AR process - ACF decays slowly
  x_ar <- numeric(100)
  x_ar[1] <- rnorm(1)
  for (i in 2:100) {
    x_ar[i] <- 0.8 * x_ar[i - 1] + rnorm(1)
  }
  result_ar <- ts_misc(x_ar)

  # AR should have later zero crossing than white noise
  expect_true(is.numeric(result_ar$firstzero_ac))
  expect_true(result_ar$firstzero_ac >= result_wn$firstzero_ac)
})

test_that("Zero proportion works correctly", {
  # Series with known zeros
  x_zeros <- c(0, 0, 1, 2, 0, 3, 0, 0, 4, 5)
  result <- ts_misc(x_zeros)

  # Should be exactly 0.5
  expect_equal(result$zero_proportion, 0.5)

  # Series with no zeros
  x_no_zeros <- c(1, 2, 3, 4, 5)
  result_no_zeros <- ts_misc(x_no_zeros)
  expect_equal(result_no_zeros$zero_proportion, 0)

  # Series with all zeros
  x_all_zeros <- rep(0, 10)
  result_all_zeros <- ts_misc(x_all_zeros)
  expect_equal(result_all_zeros$zero_proportion, 1)
})

test_that("Zero proportion with tolerance works", {
  # Series with near-zeros
  x_near <- c(1e-9, 1e-10, 1, 2, 1e-8, 3)

  # Default tolerance 1e-8
  result_default <- cpp_zero_proportion(x_near)
  expect_equal(result_default, 2 / 6)  # 1e-9 and 1e-10 are < 1e-8

  # Stricter tolerance
  result_strict <- cpp_zero_proportion(x_near, tol = 1e-9)
  expect_equal(result_strict, 1 / 6)  # Only 1e-10 is < 1e-9

  # More lenient tolerance
  result_lenient <- cpp_zero_proportion(x_near, tol = 1e-7)
  expect_equal(result_lenient, 3 / 6)  # 1e-9, 1e-10, and 1e-8 are all < 1e-7
})

test_that("Std of first derivative works", {
  set.seed(789)

  x <- rnorm(100)
  result <- ts_misc(x)

  # Compare with R's sd(diff())
  r_sd <- sd(diff(x))
  expect_equal(result$std1st_der, r_sd, tolerance = 1e-10)
})

test_that("Std of first derivative for smooth vs noisy data", {
  set.seed(111)

  # Smooth data (low derivative)
  x_smooth <- seq(0, 10, length.out = 100)
  result_smooth <- ts_misc(x_smooth)

  # Noisy data (high derivative)
  x_noisy <- rnorm(100, sd = 10)
  result_noisy <- ts_misc(x_noisy)

  # Noisy should have higher std of derivative
  expect_true(result_noisy$std1st_der > result_smooth$std1st_der)
})

test_that("Misc features handle short series", {
  # Too short for some features
  x_short <- rnorm(5)

  result <- ts_misc(x_short)

  # Should still compute
  expect_equal(length(result), 3)
  expect_true(is.numeric(result$firstzero_ac))
  expect_true(is.numeric(result$zero_proportion))
  expect_true(is.numeric(result$std1st_der))
})

test_that("Misc features integrate with ts_features", {
  set.seed(222)

  x <- rnorm(100)

  # Test with 'misc' subset
  misc_feats <- ts_features(x, features = "misc")
  expect_equal(length(misc_feats), 3)
  expect_true("firstzero_ac" %in% names(misc_feats))
  expect_true("zero_proportion" %in% names(misc_feats))
  expect_true("std1st_der" %in% names(misc_feats))

  # Test with 'all'
  all_feats <- ts_features(x, features = "all")
  expect_true("firstzero_ac" %in% names(all_feats))
  expect_true("zero_proportion" %in% names(all_feats))
  expect_true("std1st_der" %in% names(all_feats))
  expect_true(length(all_feats) >= 121)  # Should have at least 121 features now
})

test_that("First zero AC returns 0 for no crossing", {
  # Perfectly autocorrelated series
  x_perfect <- rep(1, 100)
  result <- cpp_firstzero_ac(x_perfect)

  # Should return 0 (no zero crossing)
  expect_equal(result, 0)
})

test_that("First zero AC works with different max_lag", {
  set.seed(333)

  x <- rnorm(100)

  # Small max_lag
  result_small <- cpp_firstzero_ac(x, max_lag = 10)
  expect_true(result_small >= 0)

  # Large max_lag
  result_large <- cpp_firstzero_ac(x, max_lag = 50)
  expect_true(result_large >= 0)

  # If found within small window, should match
  if (result_small > 0 && result_small <= 10) {
    expect_equal(result_small, result_large)
  }
})

test_that("Misc features handle constant series", {
  x_const <- rep(5, 100)

  result <- ts_misc(x_const)

  # Zero proportion should be 0 (no zeros)
  expect_equal(result$zero_proportion, 0)

  # Std of first derivative should be 0 or very small
  expect_true(result$std1st_der < 1e-10 || is.na(result$std1st_der))

  # First zero AC should be 0 (no crossing for constant ACF)
  expect_equal(result$firstzero_ac, 0)
})

test_that("Misc features handle NA values", {
  set.seed(444)

  # Series with NAs (will be removed by ts_features)
  x <- rnorm(100)
  x[c(10, 50, 80)] <- NA

  # ts_features removes NAs before calling ts_misc
  all_feats <- ts_features(x, features = "misc")

  # Should handle NA values
  expect_equal(length(all_feats), EXPECTED_MISC_FEATURES)
  expect_true(is.numeric(all_feats$firstzero_ac) || is.na(all_feats$firstzero_ac))
  expect_true(is.numeric(all_feats$zero_proportion) || is.na(all_feats$zero_proportion))
  expect_true(is.numeric(all_feats$std1st_der) || is.na(all_feats$std1st_der))
})

test_that("Zero proportion handles NAs correctly", {
  # Test with NAs in input
  x_with_na <- c(0, NA, 1, 0, NA, 2, 0)

  result <- cpp_zero_proportion(x_with_na)

  # Should count 3 zeros out of 5 non-NA values = 0.6
  expect_equal(result, 0.6)
})

test_that("Std first derivative matches R exactly", {
  set.seed(555)

  # Test multiple series
  for (i in 1:10) {
    x <- rnorm(100)
    our_result <- cpp_std1st_der(x)
    r_result <- sd(diff(x))

    expect_equal(our_result, r_result, tolerance = 1e-10)
  }
})

test_that("Misc features are consistent across runs", {
  set.seed(666)
  x <- rnorm(100)

  result1 <- ts_misc(x)
  result2 <- ts_misc(x)

  # Should get identical results
  expect_equal(result1$firstzero_ac, result2$firstzero_ac)
  expect_equal(result1$zero_proportion, result2$zero_proportion)
  expect_equal(result1$std1st_der, result2$std1st_der, tolerance = 1e-10)
})

test_that("First zero AC detects seasonal patterns", {
  set.seed(777)

  # Seasonal data with period 12
  t <- 1:120
  x_seasonal <- sin(2 * pi * t / 12) + rnorm(120, sd = 0.1)

  result <- ts_misc(x_seasonal)

  # ACF should exist and be numeric
  expect_true(is.numeric(result$firstzero_ac))
  expect_true(result$firstzero_ac >= 0)

  # For comparison, white noise should cross zero quickly
  x_wn <- rnorm(120)
  result_wn <- ts_misc(x_wn)

  # Seasonal pattern typically has different crossing behavior than white noise
  expect_true(is.numeric(result_wn$firstzero_ac))
})

test_that("Zero proportion detects sparse series", {
  # Sparse series (mostly zeros)
  x_sparse <- rep(0, 100)
  x_sparse[c(10, 30, 50, 70, 90)] <- rnorm(5)

  result <- ts_misc(x_sparse)

  # Should have high zero proportion
  expect_true(result$zero_proportion >= 0.9)
})

test_that("Std first derivative detects volatility", {
  set.seed(888)

  # Low volatility
  x_low_vol <- rnorm(100, sd = 0.1)
  result_low <- ts_misc(x_low_vol)

  # High volatility
  x_high_vol <- rnorm(100, sd = 10)
  result_high <- ts_misc(x_high_vol)

  # High volatility should have higher derivative std
  expect_true(result_high$std1st_der > result_low$std1st_der)
})

test_that("Feature count is correct after Phase 3 Batch 6", {
  set.seed(999)
  x <- rnorm(150)

  all_feats <- ts_features(x, features = "all")

  # Should have 121 total features (118 + 3)
  expect_equal(length(all_feats), EXPECTED_ALL_FEATURES)
})

test_that("First zero AC works for different ACF patterns", {
  set.seed(101010)

  # Quick decay (white noise)
  x_wn <- rnorm(100)
  result_wn <- cpp_firstzero_ac(x_wn)

  # Slow decay (AR)
  x_ar <- filter(rnorm(100), filter = 0.9, method = "recursive")
  x_ar <- as.numeric(x_ar)
  result_ar <- cpp_firstzero_ac(x_ar)

  # AR should have later crossing
  expect_true(result_ar >= result_wn)
})

test_that("All misc features have correct types", {
  set.seed(111111)
  x <- rnorm(100)

  result <- ts_misc(x)

  # firstzero_ac should be integer
  expect_true(is.numeric(result$firstzero_ac))
  expect_equal(result$firstzero_ac, as.integer(result$firstzero_ac))

  # zero_proportion should be in [0, 1]
  expect_true(result$zero_proportion >= 0)
  expect_true(result$zero_proportion <= 1)

  # std1st_der should be non-negative
  expect_true(result$std1st_der >= 0)
})

test_that("Misc features work with trending data", {
  set.seed(121212)

  # Strong trend
  t <- 1:100
  x_trend <- 0.1 * t + rnorm(100, sd = 0.5)

  result <- ts_misc(x_trend)

  # Should compute all features
  expect_true(is.numeric(result$firstzero_ac))
  expect_true(is.numeric(result$zero_proportion))
  expect_true(is.numeric(result$std1st_der))
})

test_that("Zero proportion with all near-zero values", {
  # All values very close to zero
  x_near_zero <- rnorm(100, mean = 0, sd = 1e-10)

  result <- cpp_zero_proportion(x_near_zero)

  # Should be close to 1
  expect_true(result > 0.99)
})
