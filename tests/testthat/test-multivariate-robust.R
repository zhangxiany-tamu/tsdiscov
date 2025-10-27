# Tests for PR 2.2: Robust Matrix Statistics
# Verify numerical stability for near-singular and pathological cases

test_that("robust_log_det handles near-singular matrices", {
  # Near-zero eigenvalues
  eigenvalues <- c(1.0, 0.5, 1e-12, 1e-15)

  # Should not return -Inf
  log_det <- robust_log_det(eigenvalues, floor = 1e-10)
  expect_true(is.finite(log_det))
  expect_false(is.na(log_det))

  # With floor, small eigenvalues are replaced
  # log_det ≈ log(1) + log(0.5) + log(1e-10) + log(1e-10)
  expected <- log(1) + log(0.5) + 2 * log(1e-10)
  expect_equal(log_det, expected, tolerance = 1e-10)
})

test_that("robust_log_det handles identity matrix", {
  # Identity matrix has all eigenvalues = 1
  eigenvalues <- rep(1, 5)

  log_det <- robust_log_det(eigenvalues)
  expect_equal(log_det, 0, tolerance = 1e-10)
})

test_that("robust_cond_number returns NA for singular matrices", {
  # Singular matrix (has zero eigenvalue)
  eigenvalues <- c(5.0, 2.0, 1.0, 1e-12)

  cond_num <- robust_cond_number(eigenvalues, eps = 1e-10)
  expect_true(is.na(cond_num))
})

test_that("robust_cond_number works for well-conditioned matrices", {
  # Well-conditioned (eigenvalues not too different)
  eigenvalues <- c(2.0, 1.5, 1.2, 1.0)

  cond_num <- robust_cond_number(eigenvalues, eps = 1e-10)
  expect_false(is.na(cond_num))
  expect_equal(cond_num, 2.0 / 1.0)
})

test_that("robust_cond_number works for identity matrix", {
  eigenvalues <- rep(1, 5)

  cond_num <- robust_cond_number(eigenvalues)
  expect_equal(cond_num, 1.0)
})

test_that("correlation features remain finite for near-singular correlation matrices", {
  set.seed(123)
  # Create nearly perfectly correlated data
  X <- matrix(0, nrow = 3, ncol = 100)
  base_signal <- rnorm(100)
  X[1,] <- base_signal
  X[2,] <- base_signal + rnorm(100, sd = 0.0001)  # Nearly identical
  X[3,] <- base_signal + rnorm(100, sd = 0.0001)

  # Extract correlation features
  corr_feat <- ts_features_multivariate(X, features = "correlation")

  # Should not have any -Inf or NaN
  expect_false(any(is.infinite(unlist(corr_feat))))
  expect_false(any(is.nan(unlist(corr_feat))))

  # Log-determinant should be finite (not -Inf)
  expect_true(is.finite(corr_feat$corr_log_determinant))

  # Condition number might be NA (acceptable for near-singular)
  # or a large number if not quite singular
  if (!is.na(corr_feat$corr_condition_number)) {
    expect_true(corr_feat$corr_condition_number > 1)
  }

  # Should still return 15 features
  expect_length(corr_feat, 15)
})

test_that("covariance features remain finite for near-singular covariance matrices", {
  set.seed(456)
  # Create nearly collinear data
  X <- matrix(0, nrow = 3, ncol = 100)
  X[1,] <- rnorm(100, sd = 10)
  X[2,] <- X[1,] * 0.99 + rnorm(100, sd = 0.01)  # Nearly proportional
  X[3,] <- X[1,] * 0.98 + rnorm(100, sd = 0.01)

  # Extract covariance features (suppress expected warning)
  cov_feat <- suppressWarnings(
    ts_features_multivariate(X, features = "covariance")
  )

  # Should not have any -Inf or NaN
  expect_false(any(is.infinite(unlist(cov_feat))))
  expect_false(any(is.nan(unlist(cov_feat))))

  # Log-determinant should be finite
  expect_true(is.finite(cov_feat$cov_log_determinant))

  # Condition number might be NA for near-singular
  if (!is.na(cov_feat$cov_condition_number)) {
    expect_true(cov_feat$cov_condition_number > 1)
  }

  # Should return 5 features (updated from 6)
  expect_length(cov_feat, 5)
})

test_that("all multivariate features remain finite for pathological data", {
  set.seed(789)
  # Perfectly correlated data (worst case)
  X <- matrix(0, nrow = 4, ncol = 100)
  base <- rnorm(100)
  for (i in 1:4) {
    X[i,] <- base + rnorm(100, sd = 1e-6)
  }

  # Extract all features
  all_feat <- ts_features_multivariate(X, features = "all")

  # Count non-NA finite values
  non_na_vals <- unlist(all_feat)[!is.na(unlist(all_feat))]
  finite_vals <- non_na_vals[is.finite(non_na_vals)]

  # Most values should be finite (some NA is acceptable, no Inf)
  expect_true(length(finite_vals) >= 42)  # At least 42/61 finite

  # Explicitly check no -Inf
  expect_false(any(sapply(all_feat, function(x) identical(x, -Inf))))
  expect_false(any(sapply(all_feat, function(x) identical(x, Inf))))
})

test_that("identity correlation matrix produces expected results", {
  set.seed(101)
  # Independent series (correlation matrix ≈ I)
  X <- matrix(rnorm(4 * 100), nrow = 4)

  corr_feat <- ts_features_multivariate(X, features = "correlation")

  # Log-determinant should be near 0 (det(I) = 1, log(1) = 0)
  expect_true(abs(corr_feat$corr_log_determinant) < 2)

  # Condition number should be near 1 (well-conditioned)
  if (!is.na(corr_feat$corr_condition_number)) {
    expect_true(corr_feat$corr_condition_number < 5)
  }

  # Mean correlation should be near 0
  expect_true(abs(corr_feat$corr_mean) < 0.3)
})

test_that("nearly constant series produce valid features", {
  set.seed(202)
  # One nearly constant series (very low variance)
  X <- matrix(0, nrow = 3, ncol = 100)
  X[1,] <- rnorm(100, sd = 1)
  X[2,] <- rnorm(100, sd = 1)
  X[3,] <- 5 + rnorm(100, sd = 0.001)  # Nearly constant

  # Should not crash (main goal of this test)
  feat <- ts_features_multivariate(X, features = "all")

  # Should return 61 features
  expect_length(feat, 61)

  # All should be numeric scalars (NA is acceptable for some)
  expect_true(all(sapply(feat, length) == 1))
  expect_true(all(sapply(feat, function(x) is.numeric(x) || is.na(x))))
})

test_that("robust functions match expected behavior for well-conditioned cases", {
  set.seed(303)
  # Well-conditioned correlation matrix
  X <- matrix(rnorm(4 * 150), nrow = 4)
  X[2,] <- 0.5 * X[1,] + 0.5 * rnorm(150)  # Moderate correlation

  # Extract features
  corr_feat <- ts_features_multivariate(X, features = "correlation")

  # All features should be finite and non-NA
  expect_true(all(is.finite(unlist(corr_feat))))
  expect_false(any(is.na(unlist(corr_feat))))
})

test_that("very small eigenvalues are handled correctly", {
  # Simulate eigenvalues from a nearly singular matrix
  eigenvalues <- c(3.0, 1.0, 0.1, 1e-20, 1e-30)

  # Log-determinant with floor
  log_det <- robust_log_det(eigenvalues, floor = 1e-10)
  expect_true(is.finite(log_det))

  # Should approximately equal: log(3) + log(1) + log(0.1) + 2*log(1e-10)
  expected <- log(3) + log(1) + log(0.1) + 2 * log(1e-10)
  expect_equal(log_det, expected, tolerance = 1e-9)

  # Condition number should be NA (smallest eigenvalue too small)
  cond_num <- robust_cond_number(eigenvalues, eps = 1e-10)
  expect_true(is.na(cond_num))
})

test_that("total feature count is correct after robust changes", {
  set.seed(404)
  X <- matrix(rnorm(4 * 100), nrow = 4)

  # Extract all features
  all_feat <- ts_features_multivariate(X, features = "all")

  # Should return 61 features total:
  # PCA: 15 + Correlation: 15 + Covariance: 5 + Sync: 8 + Diversity: 7 + Total Correlation: 2 = 61
  expect_length(all_feat, 61)
})
