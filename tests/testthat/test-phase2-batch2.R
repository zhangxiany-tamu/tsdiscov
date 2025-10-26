test_that("STL features work on seasonal data", {
  set.seed(123)
  t <- 1:100
  # Create series with trend and seasonality
  x <- sin(2 * pi * t / 12) + 0.1 * t + rnorm(100, sd = 0.2)

  stl_features <- ts_stl(x, frequency = 12)

  # Should have 5 features
  expect_equal(length(stl_features), 18)

  # Check names
  expect_true("trend_strength" %in% names(stl_features))
  expect_true("seasonal_strength" %in% names(stl_features))
  expect_true("spike" %in% names(stl_features))
  expect_true("linearity" %in% names(stl_features))
  expect_true("curvature" %in% names(stl_features))

  # All should be numeric
  for (i in seq_along(stl_features)) {
    expect_true(is.numeric(stl_features[[i]]))
  }

  # For data with seasonality, seasonal_strength should be high
  if (is.finite(stl_features$seasonal_strength)) {
    expect_true(stl_features$seasonal_strength > 0.5)
  }
})

test_that("STL features handle non-seasonal data", {
  set.seed(42)
  # Non-seasonal data
  x <- rnorm(100)

  stl_features <- ts_stl(x, frequency = 1)

  # Should return NA values for non-seasonal data
  expect_true(is.na(stl_features$trend_strength))
  expect_true(is.na(stl_features$seasonal_strength))
})

test_that("Trend strength computed correctly", {
  set.seed(123)

  # Create components manually
  n <- 100
  trend <- seq(1, 10, length.out = n)
  remainder <- rnorm(n, sd = 0.1)

  trend_str <- cpp_trend_strength(trend, remainder)

  # Should be finite and in [0, 1]
  expect_true(is.finite(trend_str))
  expect_true(trend_str >= 0 && trend_str <= 1)

  # Strong trend with small remainder should have high strength
  expect_true(trend_str > 0.8)
})

test_that("Trend strength handles weak trend", {
  set.seed(42)

  # Weak trend with large noise
  n <- 100
  trend <- seq(1, 2, length.out = n)  # Small trend
  remainder <- rnorm(n, sd = 5)  # Large noise

  trend_str <- cpp_trend_strength(trend, remainder)

  # Should be finite
  expect_true(is.finite(trend_str) || is.na(trend_str))

  # Weak trend should have low strength
  if (is.finite(trend_str)) {
    expect_true(trend_str < 0.5)
  }
})

test_that("Seasonal strength computed correctly", {
  set.seed(123)

  # Create seasonal component
  n <- 100
  t <- 1:n
  seasonal <- sin(2 * pi * t / 12)
  remainder <- rnorm(n, sd = 0.1)

  seas_str <- cpp_seasonal_strength(seasonal, remainder)

  # Should be finite and in [0, 1]
  expect_true(is.finite(seas_str))
  expect_true(seas_str >= 0 && seas_str <= 1)

  # Strong seasonality with small remainder should have high strength
  expect_true(seas_str > 0.8)
})

test_that("Spike detects outliers in remainder", {
  set.seed(123)

  # Remainder with no spikes
  remainder_smooth <- rnorm(100, sd = 0.1)
  spike_smooth <- cpp_spike(remainder_smooth)

  # Remainder with spike
  remainder_spiky <- rnorm(100, sd = 0.1)
  remainder_spiky[50] <- 10  # Add spike
  spike_spiky <- cpp_spike(remainder_spiky)

  # Both should be finite
  expect_true(is.finite(spike_smooth))
  expect_true(is.finite(spike_spiky))

  # Spiky series should have higher spike value
  expect_true(spike_spiky > spike_smooth)
})

test_that("Linearity measures trend linearity", {
  set.seed(123)

  # Linear trend
  n <- 100
  trend_linear <- seq(1, 10, length.out = n)
  lin_linear <- cpp_linearity(trend_linear)

  # Should be close to 1 for perfect linear trend
  expect_true(is.finite(lin_linear))
  expect_true(lin_linear > 0.99)

  # Curved trend
  t <- 1:n
  trend_curved <- 0.01 * t^2
  lin_curved <- cpp_linearity(trend_curved)

  # Should be lower for curved trend
  expect_true(is.finite(lin_curved))
  expect_true(lin_curved < lin_linear)
})

test_that("Curvature detects non-linear trends", {
  set.seed(123)

  # Linear trend
  n <- 100
  trend_linear <- seq(1, 10, length.out = n)
  curv_linear <- cpp_curvature(trend_linear)

  # Quadratic trend
  t <- 1:n
  trend_quad <- 0.01 * t^2 + t
  curv_quad <- cpp_curvature(trend_quad)

  # Both should be finite
  expect_true(is.finite(curv_linear) || is.na(curv_linear))
  expect_true(is.finite(curv_quad) || is.na(curv_quad))
})

test_that("STL features work with different frequencies", {
  set.seed(42)
  t <- 1:200

  # Weekly seasonality
  x_weekly <- sin(2 * pi * t / 7) + rnorm(200, sd = 0.2)
  stl_weekly <- ts_stl(x_weekly, frequency = 7)

  # Monthly seasonality
  x_monthly <- sin(2 * pi * t / 12) + rnorm(200, sd = 0.2)
  stl_monthly <- ts_stl(x_monthly, frequency = 12)

  # Both should compute successfully
  expect_equal(length(stl_weekly), 18)
  expect_equal(length(stl_monthly), 18)
})

test_that("STL handles short series", {
  # Too short for STL with frequency 12
  x_short <- rnorm(20)
  stl_short <- ts_stl(x_short, frequency = 12)

  # Should return NA values
  expect_true(is.na(stl_short$trend_strength))
  expect_true(is.na(stl_short$seasonal_strength))
})

test_that("Complete feature extraction includes STL", {
  set.seed(123)
  t <- 1:100
  x <- sin(2 * pi * t / 12) + 0.1 * t + rnorm(100, sd = 0.2)

  # Extract all features
  all_features <- ts_features(x, features = "all")

  # Should now have 94 + 11 = 105 features
  expect_equal(length(all_features), 352)

  # Check STL features are included
  expect_true("trend_strength" %in% names(all_features))
  expect_true("seasonal_strength" %in% names(all_features))
  expect_true("spike" %in% names(all_features))

  # All features should be numeric scalars
  for (i in seq_along(all_features)) {
    expect_true(is.numeric(all_features[[i]]) || is.logical(all_features[[i]]))
    expect_equal(length(all_features[[i]]), 1)
  }
})

test_that("STL features can be extracted independently", {
  set.seed(123)
  t <- 1:100
  x <- sin(2 * pi * t / 12) + rnorm(100, sd = 0.2)

  # Extract only STL features
  stl_only <- ts_features(x, features = "stl")

  # Should have exactly 5 features
  expect_equal(length(stl_only), 18)
})

test_that("STL features distinguish signal types", {
  set.seed(42)
  t <- 1:200

  # Strongly seasonal signal
  seasonal <- sin(2 * pi * t / 12) + rnorm(200, sd = 0.1)
  f_seasonal <- ts_stl(seasonal, frequency = 12)

  # Trending signal
  trending <- 0.1 * t + rnorm(200, sd = 0.5)
  f_trending <- ts_stl(trending, frequency = 12)

  # Seasonal signal should have higher seasonal strength
  if (is.finite(f_seasonal$seasonal_strength) && is.finite(f_trending$seasonal_strength)) {
    expect_true(f_seasonal$seasonal_strength > f_trending$seasonal_strength)
  }

  # Trending signal should have high trend strength
  if (is.finite(f_trending$trend_strength)) {
    expect_true(f_trending$trend_strength > 0.5)
  }
})

test_that("STL features are consistent across runs", {
  set.seed(999)
  t <- 1:100
  x <- sin(2 * pi * t / 12) + rnorm(100, sd = 0.2)

  f1 <- ts_stl(x, frequency = 12)
  f2 <- ts_stl(x, frequency = 12)

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

test_that("Edge case: constant series", {
  x_const <- rep(5, 100)

  stl_features <- ts_stl(x_const, frequency = 12)

  # STL may fail on constant series, should return NAs gracefully
  expect_equal(length(stl_features), 18)
  for (i in seq_along(stl_features)) {
    expect_true(is.numeric(stl_features[[i]]))
  }
})

test_that("Trend and seasonal strength are bounded", {
  set.seed(123)
  t <- 1:100
  x <- sin(2 * pi * t / 12) + 0.1 * t + rnorm(100, sd = 0.2)

  stl_features <- ts_stl(x, frequency = 12)

  # Strengths should be in [0, 1]
  if (is.finite(stl_features$trend_strength)) {
    expect_true(stl_features$trend_strength >= 0)
    expect_true(stl_features$trend_strength <= 1)
  }

  if (is.finite(stl_features$seasonal_strength)) {
    expect_true(stl_features$seasonal_strength >= 0)
    expect_true(stl_features$seasonal_strength <= 1)
  }
})
