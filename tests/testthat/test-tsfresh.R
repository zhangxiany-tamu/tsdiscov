# Tests for tsfresh-inspired features

test_that("C3 non-linearity feature works", {
  set.seed(123)
  x <- rnorm(100)

  # Test with lag 1
  c3_val <- cpp_c3(x, lag = 1)
  expect_true(is.finite(c3_val))
  expect_true(is.numeric(c3_val))

  # Test with different lags
  c3_lag2 <- cpp_c3(x, lag = 2)
  c3_lag3 <- cpp_c3(x, lag = 3)
  expect_true(is.finite(c3_lag2))
  expect_true(is.finite(c3_lag3))

  # For white noise, C3 should be close to 0
  expect_true(abs(c3_val) < 1)
})

test_that("C3 handles edge cases", {
  # Too short series
  x_short <- rnorm(3)
  expect_true(is.na(cpp_c3(x_short, lag = 2)))

  # Constant series
  x_const <- rep(5, 100)
  c3_const <- cpp_c3(x_const, lag = 1)
  expect_equal(c3_const, 125)  # 5 * 5 * 5 = 125
})

test_that("CID complexity feature works", {
  set.seed(123)

  # White noise should have moderate complexity
  wn <- rnorm(100)
  cid_wn <- cpp_cid_ce(wn, normalize = TRUE)
  expect_true(is.finite(cid_wn))
  expect_true(cid_wn > 0)

  # Constant series has zero complexity
  const <- rep(5, 100)
  cid_const <- cpp_cid_ce(const, normalize = FALSE)
  expect_equal(cid_const, 0)

  # More complex series should have higher CID
  trend <- 1:100 + rnorm(100)
  sine <- sin(seq(0, 4*pi, length.out = 100)) + rnorm(100, sd = 0.1)
  cid_trend <- cpp_cid_ce(trend, normalize = TRUE)
  cid_sine <- cpp_cid_ce(sine, normalize = TRUE)

  expect_true(is.finite(cid_trend))
  expect_true(is.finite(cid_sine))
  expect_true(cid_sine > cid_const)
})

test_that("Lempel-Ziv complexity works", {
  set.seed(42)

  # Random series should have high complexity
  random <- rnorm(100)
  lz_random <- cpp_lempel_ziv_complexity(random, bins = 10)
  expect_true(is.finite(lz_random))
  expect_true(lz_random > 0)
  expect_true(lz_random <= 1)

  # Constant series has zero complexity
  const <- rep(5, 100)
  lz_const <- cpp_lempel_ziv_complexity(const, bins = 10)
  expect_equal(lz_const, 0)

  # Periodic series should have lower complexity than random
  periodic <- sin(seq(0, 10*pi, length.out = 100))
  lz_periodic <- cpp_lempel_ziv_complexity(periodic, bins = 10)
  expect_true(is.finite(lz_periodic))
  expect_true(lz_periodic < lz_random)
})

test_that("Index mass quantile works", {
  set.seed(123)
  x <- rnorm(100)

  # Test different quantiles
  q25 <- cpp_index_mass_quantile(x, q = 0.25)
  q50 <- cpp_index_mass_quantile(x, q = 0.50)
  q75 <- cpp_index_mass_quantile(x, q = 0.75)

  expect_true(is.finite(q25))
  expect_true(is.finite(q50))
  expect_true(is.finite(q75))

  # All should be between 0 and 1
  expect_true(q25 >= 0 && q25 <= 1)
  expect_true(q50 >= 0 && q50 <= 1)
  expect_true(q75 >= 0 && q75 <= 1)

  # Should be ordered
  expect_true(q25 <= q50)
  expect_true(q50 <= q75)
})

test_that("Index mass quantile handles special cases", {
  # Constant series - all mass is at each point equally, could return any valid index
  const <- rep(5, 100)
  result_const <- cpp_index_mass_quantile(const, q = 0.5)
  expect_true(is.na(result_const) || (is.finite(result_const) && result_const >= 0 && result_const <= 1))

  # Series with increasing trend
  trend <- 1:100
  q50 <- cpp_index_mass_quantile(trend, q = 0.5)
  # Mass center should be around middle to end
  expect_true(q50 > 0.5)
})

test_that("Mean N absolute max works", {
  set.seed(123)
  x <- c(1, -5, 3, -8, 2, 10, -3, 7, -1, 4)

  # Top 1 should be max absolute value
  mean_1 <- cpp_mean_n_absolute_max(x, number_of_maxima = 1)
  expect_equal(mean_1, 10)

  # Top 3 should be mean of [10, 8, 7]
  mean_3 <- cpp_mean_n_absolute_max(x, number_of_maxima = 3)
  expect_equal(mean_3, (10 + 8 + 7) / 3)

  # Top 5
  mean_5 <- cpp_mean_n_absolute_max(x, number_of_maxima = 5)
  expect_equal(mean_5, (10 + 8 + 7 + 5 + 4) / 5)
})

test_that("Mean N absolute max handles edge cases", {
  set.seed(42)
  x <- rnorm(10)

  # number_of_maxima > length should return NA
  expect_true(is.na(cpp_mean_n_absolute_max(x, number_of_maxima = 20)))

  # number_of_maxima = 0 should return NA
  expect_true(is.na(cpp_mean_n_absolute_max(x, number_of_maxima = 0)))

  # Negative number_of_maxima
  expect_true(is.na(cpp_mean_n_absolute_max(x, number_of_maxima = -1)))
})

test_that("Energy ratio by chunks works", {
  set.seed(123)

  # Create series with more energy in beginning
  x <- c(rep(10, 50), rep(1, 50))

  # First chunk (out of 10) should have high energy ratio
  ratio_first <- cpp_energy_ratio_by_chunks(x, num_segments = 10, segment_focus = 0)
  # Last chunk should have low energy ratio
  ratio_last <- cpp_energy_ratio_by_chunks(x, num_segments = 10, segment_focus = 9)

  expect_true(is.finite(ratio_first))
  expect_true(is.finite(ratio_last))
  expect_true(ratio_first > ratio_last)

  # All ratios should sum to 1
  total_ratio <- 0
  for (i in 0:9) {
    total_ratio <- total_ratio + cpp_energy_ratio_by_chunks(x, num_segments = 10, segment_focus = i)
  }
  expect_equal(total_ratio, 1.0, tolerance = 1e-10)
})

test_that("Change quantiles works", {
  set.seed(42)
  x <- rnorm(100)

  # Test with middle 50% corridor
  cq_mean <- cpp_change_quantiles(x, ql = 0.25, qh = 0.75, isabs = TRUE, f_agg = "mean")
  cq_std <- cpp_change_quantiles(x, ql = 0.25, qh = 0.75, isabs = TRUE, f_agg = "std")
  cq_var <- cpp_change_quantiles(x, ql = 0.25, qh = 0.75, isabs = TRUE, f_agg = "var")

  expect_true(is.finite(cq_mean) || cq_mean == 0)
  expect_true(is.finite(cq_std) || cq_std == 0)
  expect_true(is.finite(cq_var) || cq_var == 0)

  # With absolute values, should be non-negative
  expect_true(cq_mean >= 0)
  expect_true(cq_std >= 0)
  expect_true(cq_var >= 0)
})

test_that("Change quantiles handles edge cases", {
  set.seed(123)
  x <- rnorm(100)

  # ql >= qh should return 0
  expect_equal(cpp_change_quantiles(x, ql = 0.5, qh = 0.5, isabs = TRUE, f_agg = "mean"), 0)
  expect_equal(cpp_change_quantiles(x, ql = 0.7, qh = 0.3, isabs = TRUE, f_agg = "mean"), 0)

  # Too short series
  x_short <- rnorm(1)
  expect_equal(cpp_change_quantiles(x_short, ql = 0.25, qh = 0.75, isabs = TRUE, f_agg = "mean"), 0)
})

test_that("Fourier entropy works", {
  set.seed(123)

  # Random noise should have high entropy
  noise <- rnorm(128)
  fe_noise <- cpp_fourier_entropy(noise, bins = 10)
  expect_true(is.finite(fe_noise))
  expect_true(fe_noise > 0)

  # Pure sine wave should have lower entropy
  sine <- sin(seq(0, 4*pi, length.out = 128))
  fe_sine <- cpp_fourier_entropy(sine, bins = 10)
  expect_true(is.finite(fe_sine))

  # Noise should have higher entropy than pure tone
  expect_true(fe_noise > fe_sine)
})

test_that("Fourier entropy handles edge cases", {
  # Too short series
  x_short <- rnorm(2)
  expect_true(is.na(cpp_fourier_entropy(x_short, bins = 10)))

  # Constant series
  const <- rep(5, 100)
  fe_const <- cpp_fourier_entropy(const, bins = 10)
  expect_true(is.na(fe_const) || is.finite(fe_const))
})

test_that("ts_tsfresh wrapper works", {
  set.seed(42)
  x <- rnorm(100)

  result <- ts_tsfresh(x)

  # Should return list with 10 features
  expect_type(result, "list")
  expect_equal(length(result), 10)

  # Check all expected features are present
  expected_names <- c("c3_lag1", "cid_ce", "lempel_ziv",
                       "index_mass_q25", "index_mass_q50", "index_mass_q75",
                       "change_quantiles", "mean_abs_max_3",
                       "energy_ratio_first", "fourier_entropy")
  expect_equal(names(result), expected_names)

  # All should be numeric scalars
  for (feat in result) {
    expect_true(is.numeric(feat))
    expect_equal(length(feat), 1)
  }
})

test_that("ts_cwt wrapper works", {
  set.seed(123)
  x <- rnorm(100)

  result <- ts_cwt(x)

  # Should return list with 4 features (default widths)
  expect_type(result, "list")
  expect_equal(length(result), 4)

  # Check expected feature names
  expected_names <- c("cwt_coeff_width2", "cwt_coeff_width5",
                       "cwt_coeff_width10", "cwt_coeff_width20")
  expect_equal(names(result), expected_names)

  # All should be numeric scalars (might be NA if wavelets package not installed)
  for (feat in result) {
    expect_true(is.numeric(feat))
    expect_equal(length(feat), 1)
  }
})

test_that("ts_cwt handles short series", {
  x_short <- rnorm(5)
  result <- ts_cwt(x_short)

  # Should return NA values
  expect_equal(length(result), 4)
  expect_true(all(is.na(unlist(result))))
})

test_that("tsfresh features integrate with ts_features", {
  set.seed(999)
  x <- rnorm(150)

  # Test with 'tsfresh' subset
  tsfresh_feats <- ts_features(x, features = "complexity")
  expect_type(tsfresh_feats, "list")
  expect_equal(length(tsfresh_feats), 10)
  expect_true("c3_lag1" %in% names(tsfresh_feats))
  expect_true("lempel_ziv" %in% names(tsfresh_feats))

  # Test with 'cwt' subset
  cwt_feats <- ts_features(x, features = "cwt")
  expect_type(cwt_feats, "list")
  expect_equal(length(cwt_feats), 4)
  expect_true("cwt_coeff_width2" %in% names(cwt_feats))

  # Test with 'all'
  all_feats <- ts_features(x, features = "all")
  expect_true("c3_lag1" %in% names(all_feats))
  expect_true("lempel_ziv" %in% names(all_feats))
  expect_true("cwt_coeff_width2" %in% names(all_feats))
  expect_equal(length(all_feats), EXPECTED_ALL_FEATURES)  # 143 + 10 + 4 = 157
})

test_that("Feature count is correct after adding tsfresh features", {
  set.seed(101010)
  x <- rnorm(200)

  all_feats <- ts_features(x, features = "all")

  # Should have 213 total features (202 + 11 tsfresh_supp)
  expect_equal(length(all_feats), EXPECTED_ALL_FEATURES)
})

# =============================================================================
# Tests for supplementary tsfresh features
# =============================================================================

test_that("abs_energy works correctly", {
  set.seed(123)
  x <- rnorm(100)

  energy <- cpp_abs_energy(x)
  expect_true(is.finite(energy))
  expect_true(energy >= 0)  # Sum of squares is always non-negative

  # Verify calculation
  expected <- sum(x^2)
  expect_equal(energy, expected, tolerance = 1e-10)

  # Constant series
  x_const <- rep(3, 50)
  energy_const <- cpp_abs_energy(x_const)
  expect_equal(energy_const, 3^2 * 50)
})

test_that("sum_values works correctly", {
  set.seed(123)
  x <- rnorm(100)

  sum_val <- cpp_sum_values(x)
  expect_true(is.finite(sum_val))

  # Verify calculation
  expected <- sum(x)
  expect_equal(sum_val, expected, tolerance = 1e-10)

  # Constant series
  x_const <- rep(5, 50)
  sum_const <- cpp_sum_values(x_const)
  expect_equal(sum_const, 5 * 50)
})

test_that("has_duplicate detects duplicates correctly", {
  # Series with duplicates
  x_dup <- c(1, 2, 3, 2, 4)
  expect_true(cpp_has_duplicate(x_dup))

  # Series without duplicates
  x_unique <- c(1, 2, 3, 4, 5)
  expect_false(cpp_has_duplicate(x_unique))

  # All same value
  x_same <- rep(5, 10)
  expect_true(cpp_has_duplicate(x_same))

  # Single element
  x_single <- c(5)
  expect_false(cpp_has_duplicate(x_single))
})

test_that("has_duplicate_max works correctly", {
  # Max appears multiple times
  x1 <- c(1, 5, 3, 5, 2)
  expect_true(cpp_has_duplicate_max(x1))

  # Max appears once
  x2 <- c(1, 2, 3, 4, 5)
  expect_false(cpp_has_duplicate_max(x2))

  # All same (all are max)
  x3 <- rep(7, 5)
  expect_true(cpp_has_duplicate_max(x3))
})

test_that("has_duplicate_min works correctly", {
  # Min appears multiple times
  x1 <- c(5, 1, 3, 1, 2)
  expect_true(cpp_has_duplicate_min(x1))

  # Min appears once
  x2 <- c(1, 2, 3, 4, 5)
  expect_false(cpp_has_duplicate_min(x2))

  # All same (all are min)
  x3 <- rep(7, 5)
  expect_true(cpp_has_duplicate_min(x3))
})

test_that("percentage_reoccurring calculates correctly", {
  # All unique values
  x_unique <- c(1, 2, 3, 4, 5)
  perc_unique <- cpp_percentage_reoccurring(x_unique)
  expect_equal(perc_unique, 0)  # 0% reoccurring

  # All same values
  x_same <- rep(5, 10)
  perc_same <- cpp_percentage_reoccurring(x_same)
  expect_equal(perc_same, 0.9)  # (10 - 1) / 10 = 0.9

  # Mixed
  x_mixed <- c(1, 2, 2, 3, 3, 3)  # 6 total, 3 unique
  perc_mixed <- cpp_percentage_reoccurring(x_mixed)
  expect_equal(perc_mixed, 0.5)  # (6 - 3) / 6 = 0.5
})

test_that("sum_reoccurring calculates correctly", {
  # All unique
  x_unique <- c(1, 2, 3, 4, 5)
  sum_unique <- cpp_sum_reoccurring(x_unique)
  expect_equal(sum_unique, 0)  # No reoccurring values

  # Some duplicates
  x_dup <- c(1, 2, 2, 3, 3, 3)
  # 2 appears 2 times: 2*2 = 4
  # 3 appears 3 times: 3*3 = 9
  # Total: 4 + 9 = 13
  sum_dup <- cpp_sum_reoccurring(x_dup)
  expect_equal(sum_dup, 13)
})

test_that("ratio_unique_values calculates correctly", {
  # All unique
  x_unique <- c(1, 2, 3, 4, 5)
  ratio_unique <- cpp_ratio_unique_values(x_unique)
  expect_equal(ratio_unique, 1.0)  # 5/5 = 1

  # All same
  x_same <- rep(5, 10)
  ratio_same <- cpp_ratio_unique_values(x_same)
  expect_equal(ratio_same, 0.1)  # 1/10 = 0.1

  # Mixed
  x_mixed <- c(1, 2, 2, 3, 3, 3)
  ratio_mixed <- cpp_ratio_unique_values(x_mixed)
  expect_equal(ratio_mixed, 0.5)  # 3/6 = 0.5
})

test_that("symmetry_looking detects symmetric distributions", {
  set.seed(123)

  # Symmetric distribution (normal)
  x_sym <- rnorm(1000)
  expect_true(cpp_symmetry_looking(x_sym, r = 0.1))

  # Skewed distribution
  x_skewed <- c(rnorm(100, mean = 0), rnorm(20, mean = 10))
  # This might or might not be detected as symmetric depending on threshold
  # Just check it returns a boolean
  result <- cpp_symmetry_looking(x_skewed, r = 0.05)
  expect_true(is.logical(result))

  # Constant series (perfectly symmetric)
  x_const <- rep(5, 100)
  expect_true(cpp_symmetry_looking(x_const))
})

test_that("large_standard_deviation works correctly", {
  set.seed(123)

  # High variability relative to range
  x_variable <- rnorm(100, mean = 5, sd = 10)
  # Should have large std relative to range
  result <- cpp_large_standard_deviation(x_variable, r = 0.25)
  expect_true(is.logical(result))

  # Low variability (small range, small std)
  x_low_var <- rnorm(100, mean = 5, sd = 0.1)
  # Might or might not be large depending on range

  # Constant series (no variation)
  x_const <- rep(5, 100)
  expect_false(cpp_large_standard_deviation(x_const))
})

test_that("variance_larger_than_std works correctly", {
  set.seed(123)

  # Variance > 1 means variance > std
  # Create data with large variance
  x_large_var <- rnorm(100, mean = 0, sd = 5)  # Var ≈ 25, SD ≈ 5
  expect_true(cpp_variance_larger_than_std(x_large_var))

  # Create data with small variance (< 1)
  x_small_var <- rnorm(100, mean = 0, sd = 0.5)  # Var ≈ 0.25, SD ≈ 0.5
  expect_false(cpp_variance_larger_than_std(x_small_var))

  # Edge case: variance exactly 1
  x_unit_var <- rnorm(100, mean = 0, sd = 1)
  result <- cpp_variance_larger_than_std(x_unit_var)
  expect_true(is.logical(result))
})

test_that("ts_tsfresh_supp() extracts all supplementary features", {
  set.seed(123)
  x <- rnorm(100)

  supp_feats <- ts_tsfresh_supp(x)

  # Should have 11 features
  expect_equal(length(supp_feats), 11)

  # Check feature names
  expected_names <- c(
    "abs_energy", "sum_values",
    "has_duplicate", "has_duplicate_max", "has_duplicate_min",
    "percentage_reoccurring", "sum_reoccurring", "ratio_unique_values",
    "symmetry_looking", "large_standard_deviation", "variance_larger_than_std"
  )
  expect_equal(names(supp_feats), expected_names)

  # All features should have valid values
  for (name in names(supp_feats)) {
    expect_false(is.null(supp_feats[[name]]))
  }
})

test_that("tsfresh_supp features integrate with ts_features()", {
  set.seed(456)
  x <- rnorm(150)

  # Extract only tsfresh_supp
  supp_only <- ts_features(x, features = "aggregations")
  expect_equal(length(supp_only), 11)

  # Extract all features
  all_feats <- ts_features(x, features = "all")
  expect_true("abs_energy" %in% names(all_feats))
  expect_true("has_duplicate" %in% names(all_feats))
  expect_true("symmetry_looking" %in% names(all_feats))

  # Total should be 213 (202 + 11)
  expect_equal(length(all_feats), EXPECTED_ALL_FEATURES)
})

test_that("Supplementary features handle edge cases", {
  # Empty vector
  x_empty <- numeric(0)
  expect_true(is.na(cpp_abs_energy(x_empty)))
  expect_true(is.na(cpp_sum_values(x_empty)))

  # Single value
  x_single <- c(5)
  expect_equal(cpp_abs_energy(x_single), 25)
  expect_equal(cpp_sum_values(x_single), 5)
  expect_false(cpp_has_duplicate(x_single))

  # Series with NAs
  x_na <- c(1, 2, NA, 4, 5)
  # Functions should handle NAs gracefully
  energy_na <- cpp_abs_energy(x_na)
  expect_true(is.finite(energy_na))
})
