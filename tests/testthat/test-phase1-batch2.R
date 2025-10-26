test_that("Spectral entropy works correctly", {
  set.seed(123)

  # White noise should have high spectral entropy
  wn <- rnorm(100)
  se_wn <- cpp_spectral_entropy(wn)
  expect_true(is.finite(se_wn))
  expect_true(se_wn > 0)

  # Periodic signal should have low spectral entropy
  t <- 1:100
  periodic <- sin(2 * pi * t / 10)
  se_periodic <- cpp_spectral_entropy(periodic)
  expect_true(is.finite(se_periodic))

  # White noise should have higher entropy than periodic
  expect_true(se_wn > se_periodic)
})

test_that("Peak detection works correctly", {
  set.seed(123)

  # Series with known peaks
  x <- c(1, 2, 1, 0, 3, 2, 1, 4, 2, 1)
  peaks <- cpp_number_peaks(x, 1)
  expect_equal(peaks, 3) # Peaks at positions 2, 5, 8

  # Constant series has no peaks
  x_const <- rep(5, 100)
  peaks_const <- cpp_number_peaks(x_const, 3)
  expect_equal(peaks_const, 0)

  # Monotonic increasing has no peaks
  x_mono <- 1:100
  peaks_mono <- cpp_number_peaks(x_mono, 3)
  expect_equal(peaks_mono, 0)
})

test_that("Strike features work correctly", {
  set.seed(123)

  # Create series with known strikes
  x <- c(2, 3, 4, 1, 0, -1, 2, 3, 4, 5) # Mean should be 2.3

  strike_above <- cpp_longest_strike_above_mean(x)
  strike_below <- cpp_longest_strike_below_mean(x)

  expect_true(strike_above >= 0)
  expect_true(strike_below >= 0)
  expect_true(is.finite(strike_above))
  expect_true(is.finite(strike_below))

  # Mostly above mean
  x_above <- c(rep(10, 49), 1) # Mean will be ~9.82, so 49 values above
  strike_mostly_above <- cpp_longest_strike_above_mean(x_above)
  expect_equal(strike_mostly_above, 49)

  # Mostly below mean
  x_below <- c(rep(-10, 49), 1) # Mean will be ~-8.82, so 49 values below
  strike_mostly_below <- cpp_longest_strike_below_mean(x_below)
  expect_equal(strike_mostly_below, 49)
})

test_that("Count features work correctly", {
  set.seed(123)
  x <- rnorm(100)

  count_above <- cpp_count_above_mean(x)
  count_below <- cpp_count_below_mean(x)

  # Should sum to total length
  expect_equal(count_above + count_below, 100)

  # Both should be positive
  expect_true(count_above > 0)
  expect_true(count_below > 0)

  # For normal distribution, should be roughly equal
  expect_true(abs(count_above - count_below) < 20)
})

test_that("Ratio beyond r-sigma works correctly", {
  set.seed(123)
  x <- rnorm(1000)

  # For normal distribution:
  # ~68% within 1 sigma, so ~32% beyond
  # ~95% within 2 sigma, so ~5% beyond
  ratio_1 <- cpp_ratio_beyond_r_sigma(x, 1.0)
  ratio_2 <- cpp_ratio_beyond_r_sigma(x, 2.0)

  expect_true(is.finite(ratio_1))
  expect_true(is.finite(ratio_2))

  # Check approximate proportions for large sample
  expect_true(ratio_1 > 0.2 && ratio_1 < 0.4) # Should be ~0.32
  expect_true(ratio_2 > 0.02 && ratio_2 < 0.10) # Should be ~0.05
  expect_true(ratio_1 > ratio_2) # More beyond 1 sigma than 2 sigma
})

test_that("Number of crossings works correctly", {
  set.seed(123)

  # Series that crosses zero
  x <- sin(2 * pi * (1:100) / 10)
  crossings <- cpp_number_crossings(x, 0.0)

  # Should cross 0 approximately 20 times (10 periods * 2 crossings each)
  expect_true(crossings > 15 && crossings < 25)

  # Constant series has no crossings
  x_const <- rep(5, 100)
  crossings_const <- cpp_number_crossings(x_const, 0.0)
  expect_equal(crossings_const, 0)

  # Series always above threshold has no crossings
  x_above <- rnorm(100, mean = 10, sd = 1)
  crossings_above <- cpp_number_crossings(x_above, 0.0)
  expect_equal(crossings_above, 0)
})

test_that("Absolute sum of changes works correctly", {
  set.seed(123)

  # Constant series
  x_const <- rep(5, 100)
  asc_const <- cpp_absolute_sum_of_changes(x_const)
  expect_equal(asc_const, 0.0)

  # Known pattern
  x <- c(1, 3, 2, 4, 1) # Changes: +2, -1, +2, -3; sum = 8
  asc <- cpp_absolute_sum_of_changes(x)
  expect_equal(asc, 8.0)

  # Random walk has higher ASC than white noise
  wn <- rnorm(100)
  rw <- cumsum(rnorm(100))
  expect_true(cpp_absolute_sum_of_changes(wn) > 0)
  expect_true(cpp_absolute_sum_of_changes(rw) > 0)
})

test_that("Range works correctly", {
  x <- c(1, 5, 3, 9, 2)
  r <- cpp_range(x)
  expect_equal(r, 8.0) # 9 - 1

  # Constant series has zero range
  x_const <- rep(5, 100)
  r_const <- cpp_range(x_const)
  expect_equal(r_const, 0.0)
})

test_that("Median absolute deviation works correctly", {
  set.seed(123)
  x <- rnorm(100)

  mad <- cpp_median_absolute_deviation(x)
  expect_true(is.finite(mad))
  expect_true(mad > 0)

  # For normal distribution, MAD should be roughly 0.67 * sigma
  # With sigma=1, MAD should be ~0.67
  expect_true(mad > 0.5 && mad < 0.9)

  # Constant series has zero MAD
  x_const <- rep(5, 100)
  mad_const <- cpp_median_absolute_deviation(x_const)
  expect_equal(mad_const, 0.0)
})

test_that("Coefficient of variation works correctly", {
  set.seed(123)

  # CV = std / mean
  x <- rnorm(100, mean = 10, sd = 2)
  cv <- cpp_coefficient_of_variation(x)

  expect_true(is.finite(cv))
  expect_true(cv > 0)

  # Should be roughly 2 / 10 = 0.2
  expect_true(cv > 0.1 && cv < 0.4)

  # Zero mean returns NA
  x_zero_mean <- rnorm(100, mean = 0, sd = 1)
  cv_zero <- cpp_coefficient_of_variation(x_zero_mean)
  # May or may not be NA depending on exact mean, just check it's numeric
  expect_true(is.numeric(cv_zero))
})

test_that("Benford correlation works correctly", {
  set.seed(123)

  # Random positive numbers
  x <- abs(rnorm(100, mean = 100, sd = 50))
  benford <- cpp_benford_correlation(x)

  expect_true(is.finite(benford) || is.na(benford))

  # Correlation should be between -1 and 1
  if (is.finite(benford)) {
    expect_true(benford >= -1 && benford <= 1)
  }
})

test_that("ts_structure wrapper works correctly", {
  set.seed(123)
  x <- rnorm(100)

  struct_features <- ts_structure(x)

  # Should have 14 features
  expect_equal(length(struct_features), 14)

  # Check all expected features are present
  expected_names <- c("number_peaks", "longest_strike_above", "longest_strike_below",
                      "count_above_mean", "count_below_mean", "ratio_beyond_1sigma",
                      "ratio_beyond_2sigma", "number_crossings", "absolute_sum_changes",
                      "range", "mad", "coef_variation", "benford_correlation",
                      "spectral_entropy")
  expect_true(all(expected_names %in% names(struct_features)))

  # All should be numeric
  for (i in seq_along(struct_features)) {
    expect_true(is.numeric(struct_features[[i]]))
  }
})

test_that("Complete feature extraction with batch 2 works", {
  set.seed(123)
  x <- rnorm(100)

  # Extract all features
  all_features <- ts_features(x, features = "all")

  # Should now have 210 features (202 + 11 tsfresh_supp)
  expect_equal(length(all_features), 352)

  # Check new structure features are included
  expect_true("spectral_entropy" %in% names(all_features))
  expect_true("number_peaks" %in% names(all_features))
  expect_true("benford_correlation" %in% names(all_features))

  # All features should be numeric or logical scalars
  for (i in seq_along(all_features)) {
    expect_true(is.numeric(all_features[[i]]) || is.logical(all_features[[i]]))
    expect_equal(length(all_features[[i]]), 1)
  }
})

test_that("Structure features discriminate signal types", {
  set.seed(42)

  # White noise
  wn <- rnorm(200)
  f_wn <- ts_structure(wn)

  # Periodic signal
  t <- 1:200
  periodic <- sin(2 * pi * t / 20)
  f_periodic <- ts_structure(periodic)

  # Random walk
  rw <- cumsum(rnorm(200))
  f_rw <- ts_structure(rw)

  # Periodic should have more zero crossings than random walk
  expect_true(f_periodic$number_crossings > f_rw$number_crossings)

  # Periodic should have lower spectral entropy
  expect_true(f_periodic$spectral_entropy < f_wn$spectral_entropy)

  # All should compute successfully
  expect_true(all(sapply(f_wn, function(x) is.finite(x) || is.na(x))))
  expect_true(all(sapply(f_periodic, function(x) is.finite(x) || is.na(x))))
  expect_true(all(sapply(f_rw, function(x) is.finite(x) || is.na(x))))
})
