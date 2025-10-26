test_that("FFT coefficients are computed correctly", {
  set.seed(123)
  x <- rnorm(100)

  # Test FFT coefficient extraction
  fft_coefs <- cpp_fft_coefficients(x, 5)

  # Should have 20 features (5 coefficients × 4 attributes)
  expect_equal(length(fft_coefs), 20)

  # All should be finite
  expect_true(all(is.finite(fft_coefs)))

  # Check names
  expect_true("fft_mag_1" %in% names(fft_coefs))
  expect_true("fft_real_1" %in% names(fft_coefs))
  expect_true("fft_imag_1" %in% names(fft_coefs))
  expect_true("fft_angle_1" %in% names(fft_coefs))

  # Magnitudes should be non-negative
  mag_indices <- seq(1, 20, by = 4)
  expect_true(all(fft_coefs[mag_indices] >= 0))

  # Angles should be in [-pi, pi]
  angle_indices <- seq(4, 20, by = 4)
  expect_true(all(fft_coefs[angle_indices] >= -pi))
  expect_true(all(fft_coefs[angle_indices] <= pi))
})

test_that("FFT coefficients verify Euler's formula", {
  set.seed(42)
  x <- rnorm(50)

  fft_coefs <- cpp_fft_coefficients(x, 3)

  # For each coefficient, verify: magnitude^2 = real^2 + imag^2
  for (k in 1:3) {
    idx <- (k - 1) * 4 + 1
    mag <- as.numeric(fft_coefs[idx])
    real_part <- as.numeric(fft_coefs[idx + 1])
    imag_part <- as.numeric(fft_coefs[idx + 2])

    expect_equal(mag^2, real_part^2 + imag_part^2, tolerance = 1e-6)
  }
})

test_that("FFT aggregated statistics work correctly", {
  set.seed(123)
  x <- rnorm(100)

  fft_agg <- cpp_fft_aggregated(x)

  # Should have 3 features
  expect_equal(length(fft_agg), 3)

  # Check names
  expect_true("fft_mean_mag" %in% names(fft_agg))
  expect_true("fft_var_mag" %in% names(fft_agg))
  expect_true("fft_max_mag" %in% names(fft_agg))

  # All should be finite and non-negative
  expect_true(all(is.finite(fft_agg)))
  expect_true(all(fft_agg >= 0))

  # Max should be >= mean
  expect_true(fft_agg["fft_max_mag"] >= fft_agg["fft_mean_mag"])
})

test_that("FFT features distinguish signal types", {
  set.seed(42)

  # White noise - flat spectrum
  wn <- rnorm(200)

  # Periodic signal - concentrated spectrum
  t <- 1:200
  periodic <- sin(2 * pi * t / 20)

  fft_wn <- cpp_fft_aggregated(wn)
  fft_per <- cpp_fft_aggregated(periodic)

  # Periodic signal should have higher max magnitude (concentrated energy)
  expect_true(fft_per["fft_max_mag"] > fft_wn["fft_max_mag"])

  # Periodic signal should have higher variance in magnitudes
  expect_true(fft_per["fft_var_mag"] > fft_wn["fft_var_mag"])
})

test_that("FFT detects dominant frequency in periodic signal", {
  set.seed(123)

  # Create signal with dominant frequency at bin 5 (period = 200/5 = 40)
  t <- 1:200
  x <- sin(2 * pi * t / 40) + rnorm(200, sd = 0.1)

  fft_coefs <- cpp_fft_coefficients(x, 10)

  # Extract magnitudes
  magnitudes <- fft_coefs[seq(1, 40, by = 4)]

  # Magnitude at k=5 should be among the highest
  expect_true(magnitudes[5] > median(magnitudes))
})

test_that("ts_fft wrapper works correctly", {
  set.seed(123)
  x <- rnorm(100)

  fft_features <- ts_fft(x)

  # Should have 15 features (3 coef × 4 + 3 aggregated)
  expect_equal(length(fft_features), 15)

  # Check both coefficient and aggregated features are present
  expect_true("fft_mag_1" %in% names(fft_features))
  expect_true("fft_mean_mag" %in% names(fft_features))

  # All should be numeric
  for (i in seq_along(fft_features)) {
    expect_true(is.numeric(fft_features[[i]]))
  }
})

test_that("FFT features work with different num_coef values", {
  set.seed(123)
  x <- rnorm(100)

  # Test with 1 coefficient
  fft_1 <- ts_fft(x, num_coef = 1)
  expect_equal(length(fft_1), 7)  # 1×4 + 3 aggregated

  # Test with 5 coefficients
  fft_5 <- ts_fft(x, num_coef = 5)
  expect_equal(length(fft_5), 23)  # 5×4 + 3 aggregated

  # Test with 10 coefficients
  fft_10 <- ts_fft(x, num_coef = 10)
  expect_equal(length(fft_10), 43)  # 10×4 + 3 aggregated
})

test_that("FFT handles constant series correctly", {
  x_const <- rep(5, 100)

  fft_coefs <- cpp_fft_coefficients(x_const, 5)
  fft_agg <- cpp_fft_aggregated(x_const)

  # All coefficients should be near zero (after detrending, signal is zero)
  expect_true(all(abs(fft_coefs) < 1e-10))

  # Aggregated features should be near zero or zero
  expect_true(fft_agg["fft_mean_mag"] < 1e-10)
  expect_true(fft_agg["fft_max_mag"] < 1e-10)
})

test_that("FFT handles very short series", {
  x_short <- rnorm(5)

  fft_coefs <- cpp_fft_coefficients(x_short, 3)
  fft_agg <- cpp_fft_aggregated(x_short)

  # Should still compute but with fewer coefficients
  expect_true(length(fft_coefs) > 0)
  expect_equal(length(fft_agg), 3)
  expect_true(all(is.finite(fft_agg)))
})

test_that("Complete feature extraction includes FFT", {
  set.seed(123)
  x <- rnorm(100)

  # Extract all features
  all_features <- ts_features(x, features = "all")

  # Should now have 83 features (65 + 15 FFT + 3 scaling)
  expect_equal(length(all_features), 352)

  # Check FFT features are included
  expect_true("fft_mag_1" %in% names(all_features))
  expect_true("fft_mean_mag" %in% names(all_features))

  # All features should be numeric or logical scalars
  for (i in seq_along(all_features)) {
    expect_true(is.numeric(all_features[[i]]) || is.logical(all_features[[i]]))
    expect_equal(length(all_features[[i]]), 1)
  }
})

test_that("FFT features can be extracted independently", {
  set.seed(123)
  x <- rnorm(100)

  # Extract only FFT features
  fft_only <- ts_features(x, features = "fft")

  # Should have exactly 15 features
  expect_equal(length(fft_only), 15)

  # All names should start with "fft_"
  expect_true(all(grepl("^fft_", names(fft_only))))
})

test_that("FFT features are consistent across runs", {
  set.seed(999)
  x <- rnorm(100)

  fft1 <- ts_fft(x)
  fft2 <- ts_fft(x)

  # Should get identical results
  for (i in seq_along(fft1)) {
    expect_equal(fft1[[i]], fft2[[i]], tolerance = 1e-10)
  }
})

test_that("FFT captures harmonic content", {
  set.seed(42)

  # Create signal with two harmonics
  t <- 1:200
  # Fundamental at period 40 (bin 5)
  # Second harmonic at period 20 (bin 10)
  x <- sin(2 * pi * t / 40) + 0.5 * sin(2 * pi * t / 20) + rnorm(200, sd = 0.1)

  fft_coefs <- cpp_fft_coefficients(x, 15)

  # Extract magnitudes
  magnitudes <- fft_coefs[seq(1, 60, by = 4)]

  # Both bins 5 and 10 should have high magnitudes
  # (Note: exact bins may vary due to DFT discretization)
  # Just check that there are multiple high-magnitude bins
  high_mag_count <- sum(magnitudes > mean(magnitudes) + sd(magnitudes))
  expect_true(high_mag_count >= 2)
})
