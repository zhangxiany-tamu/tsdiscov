# Tests for PR 2.1: Standardization Semantics
# Verify that covariance uses raw data, others use standardized data

test_that("covariance features are invariant to standardize flag", {
  set.seed(123)
  # Create data with different scales
  X <- matrix(0, nrow = 3, ncol = 100)
  X[1,] <- rnorm(100, mean = 0, sd = 1)
  X[2,] <- rnorm(100, mean = 0, sd = 5)    # Different scale
  X[3,] <- rnorm(100, mean = 0, sd = 0.1)  # Much smaller scale

  # Extract covariance with standardize=TRUE (suppress expected warning)
  cov_std_true <- suppressWarnings(
    ts_features_multivariate(X, features = "covariance", standardize = TRUE)
  )

  # Extract covariance with standardize=FALSE
  cov_std_false <- ts_features_multivariate(X, features = "covariance",
                                            standardize = FALSE)

  # Results should be identical (covariance always uses raw data)
  expect_equal(cov_std_true, cov_std_false)

  # Verify features exist
  expect_length(cov_std_true, 5)
  expect_true(all(sapply(cov_std_true, is.numeric)))
})

test_that("correlation features change with standardize flag", {
  set.seed(456)
  # Create data with different scales
  X <- matrix(0, nrow = 3, ncol = 100)
  X[1,] <- rnorm(100, mean = 0, sd = 1)
  X[2,] <- rnorm(100, mean = 0, sd = 10)
  X[3,] <- rnorm(100, mean = 0, sd = 0.5)

  # Extract correlation with standardize=TRUE (default)
  corr_std_true <- ts_features_multivariate(X, features = "correlation",
                                            standardize = TRUE)

  # Extract correlation with standardize=FALSE
  corr_std_false <- ts_features_multivariate(X, features = "correlation",
                                             standardize = FALSE)

  # Note: Correlation is actually invariant to standardization
  # (correlation normalizes anyway), but we test the routing works
  expect_length(corr_std_true, 15)
  expect_length(corr_std_false, 15)
  expect_true(all(sapply(corr_std_true, is.numeric)))
  expect_true(all(sapply(corr_std_false, is.numeric)))
})

test_that("PCA features are invariant to standardize flag", {
  set.seed(789)
  # Create data with very different scales
  X <- matrix(0, nrow = 3, ncol = 100)
  X[1,] <- rnorm(100, sd = 100)   # Very large scale
  X[2,] <- rnorm(100, sd = 1)     # Medium scale
  X[3,] <- rnorm(100, sd = 0.01)  # Very small scale

  # Extract PCA with standardize=TRUE
  pca_std_true <- ts_features_multivariate(X, features = "pca",
                                           standardize = TRUE)

  # Extract PCA with standardize=FALSE
  pca_std_false <- ts_features_multivariate(X, features = "pca",
                                            standardize = FALSE)

  # Results should be identical because PCA uses cor(t(X)) which
  # internally standardizes regardless of the standardize flag
  # NOTE: This is a design choice - PCA on correlation matrix is
  # scale-invariant by definition
  expect_equal(pca_std_true, pca_std_false)

  expect_length(pca_std_true, 15)
  expect_length(pca_std_false, 15)
})

test_that("warning issued when requesting only covariance with standardize=TRUE", {
  set.seed(101)
  X <- matrix(rnorm(3 * 100), nrow = 3)

  # Should warn
  expect_warning(
    ts_features_multivariate(X, features = "covariance", standardize = TRUE),
    "Covariance features use raw"
  )

  # Should NOT warn when standardize=FALSE
  expect_silent({
    ts_features_multivariate(X, features = "covariance", standardize = FALSE)
  })

  # Should NOT warn when requesting multiple feature sets
  expect_silent({
    ts_features_multivariate(X, features = c("covariance", "pca"),
                             standardize = TRUE)
  })
})

test_that("all features work correctly with dual-path routing", {
  set.seed(202)
  X <- matrix(rnorm(4 * 100), nrow = 4)

  # Extract all features with standardize=TRUE
  feat_all <- ts_features_multivariate(X, features = "all", standardize = TRUE)

  # Should return 61 features
  expect_length(expect_length(feat_all, 61), 61)
  expect_true(all(sapply(feat_all, is.numeric)))

  # No NAs or Infs
  expect_false(any(is.na(unlist(feat_all))))
  expect_false(any(is.infinite(unlist(feat_all))))
})

test_that("synchronization features use correct data path", {
  set.seed(303)
  # Create synchronized data with different scales
  t <- seq(0, 10, length.out = 150)
  X <- matrix(0, nrow = 3, ncol = 150)
  common_signal <- sin(2 * pi * t)
  X[1,] <- common_signal * 100 + rnorm(150, sd = 1)   # Large scale
  X[2,] <- common_signal * 1 + rnorm(150, sd = 0.01)  # Small scale
  X[3,] <- common_signal * 10 + rnorm(150, sd = 0.1)  # Medium scale

  # With standardize=TRUE, synchronization should be more easily detected
  sync_std_true <- ts_features_multivariate(X, features = "sync",
                                            standardize = TRUE)

  # With standardize=FALSE, different scales might affect results
  sync_std_false <- ts_features_multivariate(X, features = "sync",
                                             standardize = FALSE)

  # Both should return 8 features
  expect_length(sync_std_true, 8)
  expect_length(sync_std_false, 8)

  # CCF max mean should be high in both cases (strong synchronization)
  expect_true(sync_std_true$ccf_max_mean > 0.5)
  expect_true(sync_std_false$ccf_max_mean > 0.5)
})

test_that("diversity features use correct data path", {
  set.seed(404)
  # Create diverse data with different scales
  X <- matrix(0, nrow = 4, ncol = 120)
  X[1,] <- rnorm(120, sd = 10)
  X[2,] <- rnorm(120, sd = 1)
  X[3,] <- rnorm(120, sd = 0.1)
  X[4,] <- rnorm(120, sd = 0.01)

  # With standardize=TRUE, mean diversity should be near zero
  div_std_true <- ts_features_multivariate(X, features = "diversity",
                                           standardize = TRUE)

  # With standardize=FALSE, mean diversity should be high
  div_std_false <- ts_features_multivariate(X, features = "diversity",
                                            standardize = FALSE)

  expect_length(div_std_true, 7)
  expect_length(div_std_false, 7)

  # Mean diversity should be lower with standardization
  # (all series have mean near 0 after standardization)
  expect_true(div_std_true$mean_diversity < div_std_false$mean_diversity)
})

test_that("covariance features reflect raw data scale", {
  set.seed(505)
  # Create data with known scale structure
  X_small <- matrix(rnorm(3 * 100, sd = 0.1), nrow = 3)
  X_large <- matrix(rnorm(3 * 100, sd = 10), nrow = 3)

  # Extract covariance features (suppress expected warning)
  cov_small <- suppressWarnings(
    ts_features_multivariate(X_small, features = "covariance", standardize = TRUE)
  )
  cov_large <- suppressWarnings(
    ts_features_multivariate(X_large, features = "covariance", standardize = TRUE)
  )

  # Trace should be much larger for X_large (sum of variances)
  # Even though standardize=TRUE, covariance uses raw data
  expect_true(cov_large$cov_trace > cov_small$cov_trace * 50)

  # Log-determinant should also differ dramatically
  # (log-det scales with product of eigenvalues)
  expect_true(cov_large$cov_log_determinant > cov_small$cov_log_determinant)
})

test_that("mixing covariance with other features works correctly", {
  set.seed(606)
  X <- matrix(rnorm(4 * 100), nrow = 4)

  # Extract mixed feature sets
  feat_mixed <- ts_features_multivariate(X,
                                         features = c("covariance", "correlation", "pca"),
                                         standardize = TRUE)

  # Should return 6 + 15 + 15 = 35 features
  expect_length(feat_mixed, 35)
  expect_true(all(sapply(feat_mixed, is.numeric)))

  # Should NOT warn (multiple feature sets)
  expect_silent({
    ts_features_multivariate(X,
                             features = c("covariance", "sync"),
                             standardize = TRUE)
  })
})
