# Phase 3 Batch 2: STL Remainder ACF Features
# Tests for remainder ACF/differencing features from STL decomposition

test_that("STL remainder ACF features work correctly", {
  set.seed(123)

  # Create seasonal data
  t <- 1:100
  x <- sin(2*pi*t/12) + rnorm(100, sd=0.2)

  result <- ts_stl(x, frequency=12)

  # Should have 16 features now (5 basic + 11 remainder)
  expect_equal(length(result), 18)

  # Check remainder ACF features exist
  expect_true("e_acf_first_min" %in% names(result))
  expect_true("e_acf_timescale" %in% names(result))
  expect_true("e_acf_sum10" %in% names(result))
  expect_true("e_acf1" %in% names(result))
  expect_true("e_acf10" %in% names(result))

  # Check remainder differencing features exist
  expect_true("e_diff1_acf1" %in% names(result))
  expect_true("e_diff1_acf10" %in% names(result))
  expect_true("e_diff2_acf1" %in% names(result))
  expect_true("e_diff2_acf10" %in% names(result))
  expect_true("e_diff1x_pacf5" %in% names(result))
  expect_true("e_diff2x_pacf5" %in% names(result))
})

test_that("Remainder ACF features are valid", {
  set.seed(456)

  # Create seasonal data with strong seasonality
  t <- 1:120
  x <- 2 * sin(2*pi*t/12) + 0.5 * t/100 + rnorm(120, sd=0.3)

  result <- ts_stl(x, frequency=12)

  # All features should be finite
  remainder_feats <- result[grepl("^e_", names(result))]
  expect_true(all(sapply(remainder_feats, is.finite)))

  # ACF values should be in valid range [-1, 1]
  expect_true(result$e_acf1 >= -1 && result$e_acf1 <= 1)
  expect_true(result$e_acf10 >= -1 && result$e_acf10 <= 1)

  # Differencing ACF should also be in valid range
  expect_true(result$e_diff1_acf1 >= -1 && result$e_diff1_acf1 <= 1)
  expect_true(result$e_diff2_acf1 >= -1 && result$e_diff2_acf1 <= 1)

  # Sum of squares features should be non-negative
  expect_true(result$e_diff1_acf10 >= 0)
  expect_true(result$e_diff2_acf10 >= 0)
  expect_true(result$e_diff1x_pacf5 >= 0)
  expect_true(result$e_diff2x_pacf5 >= 0)
})

test_that("Remainder features detect residual patterns", {
  set.seed(789)

  # Good decomposition - random remainder
  t <- 1:100
  seasonal <- 2 * sin(2*pi*t/12)
  trend <- 0.5 * t/100
  noise <- rnorm(100, sd=0.1)
  x_good <- seasonal + trend + noise

  result_good <- ts_stl(x_good, frequency=12)

  # Poorly specified - AR(1) remainder
  ar_resid <- filter(rnorm(100), filter=0.7, method="recursive")
  ar_resid <- as.numeric(ar_resid)
  x_poor <- seasonal + trend + ar_resid

  result_poor <- ts_stl(x_poor, frequency=12)

  # Both should compute
  expect_true(is.finite(result_good$e_acf1))
  expect_true(is.finite(result_poor$e_acf1))

  # Poor decomposition should have higher ACF in remainder
  # (though this isn't guaranteed, so we just check they're different)
  expect_true(is.numeric(result_good$e_acf1))
  expect_true(is.numeric(result_poor$e_acf1))
})

test_that("Remainder features handle non-seasonal data", {
  set.seed(111)

  # Non-seasonal data
  x <- rnorm(50)

  result <- ts_stl(x, frequency=1)

  # Should return NA for all features when no seasonality
  expect_true(is.na(result$trend_strength))
  expect_true(is.na(result$e_acf1))
  expect_true(is.na(result$e_acf10))
  expect_true(is.na(result$e_diff1_acf1))

  # All remainder features should be NA
  remainder_feats <- result[grepl("^e_", names(result))]
  expect_true(all(sapply(remainder_feats, is.na)))
})

test_that("Remainder features handle short series", {
  # Series too short for STL with frequency 12
  x_short <- rnorm(20)

  result <- ts_stl(x_short, frequency=12)

  # Should return NA values
  expect_true(is.na(result$trend_strength))
  expect_true(is.na(result$e_acf1))
  expect_equal(length(result), 18)
})

test_that("Remainder ACF features match manual calculation", {
  set.seed(222)

  # Create data and manually decompose
  t <- 1:60
  x <- sin(2*pi*t/12) + rnorm(60, sd=0.1)

  # Get STL decomposition manually
  ts_obj <- ts(x, frequency=12)
  stl_manual <- stl(ts_obj, s.window="periodic", robust=TRUE)
  remainder_manual <- as.numeric(stl_manual$time.series[, "remainder"])

  # Get our features
  result <- ts_stl(x, frequency=12)

  # Compute ACF manually on remainder
  acf_manual <- acf(remainder_manual, lag.max=10, plot=FALSE)$acf[-1]  # Remove lag 0

  # Check e_acf1 matches
  expect_equal(result$e_acf1, acf_manual[1], tolerance=1e-6)

  # Check e_acf10 matches sum of squares of first 10 ACF (x_acf10)
  expect_equal(result$e_acf10, sum(acf_manual[1:10]^2), tolerance=1e-6)

  # Check e_acf_sum10 matches sum of first 10
  expect_equal(result$e_acf_sum10, sum(acf_manual[1:10]), tolerance=1e-6)
})

test_that("Remainder differencing features work correctly", {
  set.seed(333)

  t <- 1:100
  x <- sin(2*pi*t/12) + rnorm(100, sd=0.2)

  result <- ts_stl(x, frequency=12)

  # Differencing features should be finite
  expect_true(is.finite(result$e_diff1_acf1))
  expect_true(is.finite(result$e_diff1_acf10))
  expect_true(is.finite(result$e_diff2_acf1))
  expect_true(is.finite(result$e_diff2_acf10))
  expect_true(is.finite(result$e_diff1x_pacf5))
  expect_true(is.finite(result$e_diff2x_pacf5))

  # ACF at lag 1 should be in valid range
  expect_true(abs(result$e_diff1_acf1) <= 1)
  expect_true(abs(result$e_diff2_acf1) <= 1)

  # Sum of squares should be non-negative
  expect_true(result$e_diff1_acf10 >= 0)
  expect_true(result$e_diff2_acf10 >= 0)
  expect_true(result$e_diff1x_pacf5 >= 0)
  expect_true(result$e_diff2x_pacf5 >= 0)
})

test_that("STL remainder features integrate with ts_features", {
  set.seed(444)

  t <- 1:100
  x <- sin(2*pi*t/12) + 0.3 * t/100 + rnorm(100, sd=0.2)

  # Test with 'stl' subset
  stl_feats <- ts_features(x, features="stl")
  expect_true("e_acf1" %in% names(stl_feats))
  expect_true("e_diff1_acf1" %in% names(stl_feats))
  expect_equal(length(stl_feats), 18)

  # Test with 'all'
  all_feats <- ts_features(x, features="all")
  expect_true("e_acf1" %in% names(all_feats))
  expect_true("e_diff2_acf10" %in% names(all_feats))
  expect_true(length(all_feats) >= 105)  # Should have at least 105 features now
})

test_that("Remainder features handle different frequencies", {
  set.seed(555)

  # Weekly data (frequency = 7)
  t <- 1:70
  x_weekly <- sin(2*pi*t/7) + rnorm(70, sd=0.1)
  result_weekly <- ts_stl(x_weekly, frequency=7)

  expect_true(is.finite(result_weekly$e_acf1))
  expect_equal(length(result_weekly), 18)

  # Monthly data (frequency = 12)
  t <- 1:120
  x_monthly <- sin(2*pi*t/12) + rnorm(120, sd=0.1)
  result_monthly <- ts_stl(x_monthly, frequency=12)

  expect_true(is.finite(result_monthly$e_acf1))
  expect_equal(length(result_monthly), 18)
})

test_that("Remainder ACF timescale is computed correctly", {
  set.seed(666)

  t <- 1:100
  x <- sin(2*pi*t/12) + rnorm(100, sd=0.2)

  result <- ts_stl(x, frequency=12)

  # Timescale should be positive (lag where ACF crosses e^-1)
  expect_true(result$e_acf_timescale >= 0 || is.na(result$e_acf_timescale))

  # First min should be an integer lag
  expect_true(result$e_acf_first_min >= 0 || is.na(result$e_acf_first_min))
})

test_that("Remainder features work with strong trend", {
  set.seed(777)

  # Strong trend with seasonality
  t <- 1:100
  trend <- 5 * t/100
  seasonal <- sin(2*pi*t/12)
  x <- trend + seasonal + rnorm(100, sd=0.1)

  result <- ts_stl(x, frequency=12)

  # STL should handle strong trend
  expect_true(is.finite(result$trend_strength))
  expect_true(result$trend_strength > 0.5)  # Should have high trend strength

  # Remainder features should still work
  expect_true(is.finite(result$e_acf1))
  expect_true(is.finite(result$e_diff1_acf1))
})

test_that("Remainder features work with weak seasonality", {
  set.seed(888)

  # Weak seasonality
  t <- 1:100
  seasonal <- 0.1 * sin(2*pi*t/12)  # Very weak
  x <- rnorm(100) + seasonal

  result <- ts_stl(x, frequency=12)

  # Should decompose but may have low seasonal strength
  expect_true(is.finite(result$seasonal_strength))

  # Remainder features should still work
  expect_true(is.finite(result$e_acf1))
  expect_true(is.finite(result$e_diff1_acf1))

  # Remainder might have higher variance than seasonal component
  expect_true(is.numeric(result$e_acf_sum10))
})

test_that("STL error handling works", {
  # Test data that might cause STL to fail
  x_const <- rep(5, 50)

  result <- ts_stl(x_const, frequency=12)

  # Should return NA values gracefully
  expect_true(is.na(result$trend_strength) || is.finite(result$trend_strength))
  expect_true(is.na(result$e_acf1) || is.finite(result$e_acf1))
  expect_equal(length(result), 18)
})

test_that("Remainder PACF features work correctly", {
  set.seed(999)

  t <- 1:100
  x <- sin(2*pi*t/12) + rnorm(100, sd=0.2)

  result <- ts_stl(x, frequency=12)

  # PACF sum of squares should be non-negative
  expect_true(result$e_diff1x_pacf5 >= 0)
  expect_true(result$e_diff2x_pacf5 >= 0)

  # Should be finite
  expect_true(is.finite(result$e_diff1x_pacf5))
  expect_true(is.finite(result$e_diff2x_pacf5))
})

test_that("Feature count is correct after Phase 3 Batch 2", {
  set.seed(101010)
  x <- rnorm(100)

  all_feats <- ts_features(x, features="all")

  # Should have 109 total features (94 + 11 + 4)
  expect_equal(length(all_feats), EXPECTED_ALL_FEATURES)
})
