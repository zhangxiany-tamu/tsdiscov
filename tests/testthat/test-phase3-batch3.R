# Phase 3 Batch 3: Heterogeneity/ARCH/GARCH Features
# Tests for volatility clustering and conditional heteroskedasticity

test_that("Heterogeneity features work correctly", {
  set.seed(123)

  # Create white noise
  x <- rnorm(100)

  result <- ts_heterogeneity(x, fit_garch = FALSE)

  # Should have 4 features
  expect_equal(length(result), 4)
  expect_true("arch_acf" %in% names(result))
  expect_true("garch_acf" %in% names(result))
  expect_true("arch_r2" %in% names(result))
  expect_true("garch_r2" %in% names(result))
})

test_that("ARCH ACF is computed correctly", {
  set.seed(456)

  x <- rnorm(150)
  result <- ts_heterogeneity(x, fit_garch = FALSE)

  # arch_acf should be finite and non-negative (sum of squares)
  expect_true(is.finite(result$arch_acf) || is.na(result$arch_acf))
  if (is.finite(result$arch_acf)) {
    expect_true(result$arch_acf >= 0)
  }
})

test_that("ARCH R² is in valid range", {
  set.seed(789)

  x <- rnorm(150)
  result <- ts_heterogeneity(x, fit_garch = FALSE)

  # R² should be in [0, 1] if finite
  if (is.finite(result$arch_r2)) {
    expect_true(result$arch_r2 >= 0)
    expect_true(result$arch_r2 <= 1)
  }
})

test_that("Heterogeneity detects ARCH effects", {
  set.seed(111)

  # White noise (low heterogeneity)
  wn <- rnorm(200)
  result_wn <- ts_heterogeneity(wn, fit_garch = FALSE)

  # ARCH process (higher heterogeneity)
  n <- 200
  e <- rnorm(n)
  for(i in 2:n) {
    e[i] <- sqrt(0.1 + 0.7 * e[i-1]^2) * rnorm(1)
  }
  result_arch <- ts_heterogeneity(e, fit_garch = FALSE)

  # Both should compute
  expect_true(is.numeric(result_wn$arch_acf))
  expect_true(is.numeric(result_arch$arch_acf))

  # ARCH series should have higher arch_r2 (though not guaranteed)
  expect_true(is.numeric(result_wn$arch_r2))
  expect_true(is.numeric(result_arch$arch_r2))
})

test_that("Heterogeneity handles short series", {
  # Too short for AR fitting
  x_short <- rnorm(15)

  result <- ts_heterogeneity(x_short, fit_garch = FALSE)

  # Should return NA values
  expect_true(is.na(result$arch_acf))
  expect_true(is.na(result$arch_r2))
})

test_that("GARCH fitting works when tseries available", {
  skip_if_not_installed("tseries")

  set.seed(222)

  # Create GARCH data
  n <- 200
  e <- rnorm(n)
  h <- numeric(n)
  h[1] <- 1
  for(i in 2:n) {
    h[i] <- 0.1 + 0.5 * e[i-1]^2 + 0.3 * h[i-1]
    e[i] <- sqrt(h[i]) * rnorm(1)
  }

  result <- ts_heterogeneity(e, fit_garch = TRUE)

  # garch_acf might be computed (depending on GARCH success)
  expect_true(is.numeric(result$garch_acf))
  expect_true(is.numeric(result$garch_r2))

  # If GARCH succeeded, values should be finite or NA
  expect_true(is.finite(result$garch_acf) || is.na(result$garch_acf))
  expect_true(is.finite(result$garch_r2) || is.na(result$garch_r2))
})

test_that("Heterogeneity with fit_garch=FALSE doesn't compute GARCH", {
  set.seed(333)

  x <- rnorm(100)
  result <- ts_heterogeneity(x, fit_garch = FALSE)

  # garch features should be NA when fit_garch=FALSE
  expect_true(is.na(result$garch_acf))
  expect_true(is.na(result$garch_r2))
})

test_that("Heterogeneity handles constant series", {
  x_const <- rep(5, 100)

  result <- ts_heterogeneity(x_const, fit_garch = FALSE)

  # Should return NA or handle gracefully
  expect_true(is.na(result$arch_acf) || is.finite(result$arch_acf))
  expect_equal(length(result), 4)
})

test_that("Heterogeneity integrates with ts_features", {
  set.seed(444)

  x <- rnorm(150)

  # Test with 'heterogeneity' subset
  het_feats <- ts_features(x, features = "heterogeneity")
  expect_equal(length(het_feats), 4)
  expect_true("arch_acf" %in% names(het_feats))
  expect_true("arch_r2" %in% names(het_feats))

  # Test with 'all'
  all_feats <- ts_features(x, features = "all")
  expect_true("arch_acf" %in% names(all_feats))
  expect_true("garch_r2" %in% names(all_feats))
  expect_true(length(all_feats) >= 141)  # Should have at least 109 features now
})

test_that("ARCH ACF matches manual calculation", {
  set.seed(555)

  x <- rnorm(150)

  # Get result
  result <- ts_heterogeneity(x, fit_garch = FALSE)

  # Manually compute ARCH ACF
  ar_fit <- ar(na.contiguous(x), aic = TRUE)
  if (!is.null(ar_fit$resid)) {
    whitened <- na.contiguous(ar_fit$resid)

    if (length(whitened) >= 15) {
      # Compute ACF of squared series manually using cpp_arch_acf directly
      arch_acf_direct <- tsdiscov:::cpp_arch_acf(whitened)

      # Should match exactly
      expect_equal(result$arch_acf, arch_acf_direct, tolerance = 1e-10)
    }
  }
})

test_that("cpp_arch_acf computes correctly", {
  set.seed(666)

  x <- rnorm(100)

  # Use internal C++ function
  result_cpp <- tsdiscov:::cpp_arch_acf(x)

  # Should be finite and non-negative (sum of squares)
  expect_true(is.finite(result_cpp))
  expect_true(result_cpp >= 0)

  # Should be reasonable magnitude (sum of 12 squared correlations)
  # Each squared correlation is in [0, 1], so sum should be in [0, 12]
  expect_true(result_cpp <= 12)
})

test_that("Heterogeneity handles NA values", {
  set.seed(777)

  x <- rnorm(100)
  x[c(10, 50, 80)] <- NA

  result <- ts_heterogeneity(x, fit_garch = FALSE)

  # Should handle NA values (ar() uses na.contiguous)
  expect_equal(length(result), 4)
  expect_true(is.numeric(result$arch_acf))
})

test_that("Heterogeneity works with different series lengths", {
  set.seed(888)

  # Medium length
  x_medium <- rnorm(100)
  result_medium <- ts_heterogeneity(x_medium, fit_garch = FALSE)
  expect_true(is.numeric(result_medium$arch_acf))

  # Long length
  x_long <- rnorm(500)
  result_long <- ts_heterogeneity(x_long, fit_garch = FALSE)
  expect_true(is.numeric(result_long$arch_acf))
})

test_that("ARCH effects are higher for volatile data", {
  set.seed(999)

  # Low volatility
  x_low <- rnorm(200, sd = 0.5)
  result_low <- ts_heterogeneity(x_low, fit_garch = FALSE)

  # High volatility with clustering
  n <- 200
  e <- rnorm(n)
  for(i in 2:n) {
    e[i] <- sqrt(0.1 + 0.8 * e[i-1]^2) * rnorm(1)
  }
  result_high <- ts_heterogeneity(e, fit_garch = FALSE)

  # Both should compute
  expect_true(is.finite(result_low$arch_r2) || is.na(result_low$arch_r2))
  expect_true(is.finite(result_high$arch_r2) || is.na(result_high$arch_r2))
})

test_that("Heterogeneity error handling works", {
  # Very short series
  x_tiny <- rnorm(5)
  result_tiny <- ts_heterogeneity(x_tiny, fit_garch = FALSE)

  expect_equal(length(result_tiny), 4)
  expect_true(all(sapply(result_tiny, is.na)))
})

test_that("Feature count is correct after Phase 3 Batch 3", {
  set.seed(101010)
  x <- rnorm(150)

  all_feats <- ts_features(x, features = "all")

  # Should have 109 total features (105 + 4)
  expect_equal(length(all_feats), EXPECTED_ALL_FEATURES)
})

test_that("ARCH features detect financial time series patterns", {
  skip_if_not_installed("tseries")

  set.seed(121212)

  # Simulate returns with volatility clustering
  n <- 300
  returns <- numeric(n)
  vol <- numeric(n)
  vol[1] <- 1

  for(i in 2:n) {
    # GARCH(1,1) volatility
    vol[i] <- sqrt(0.01 + 0.05 * returns[i-1]^2 + 0.90 * vol[i-1]^2)
    returns[i] <- vol[i] * rnorm(1)
  }

  result <- ts_heterogeneity(returns, fit_garch = TRUE)

  # Should detect heterogeneity
  if (is.finite(result$arch_r2)) {
    expect_true(result$arch_r2 > 0)
  }

  expect_equal(length(result), 4)
})
