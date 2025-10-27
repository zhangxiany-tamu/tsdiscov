context("Multivariate Features - Phase 1")

# Test data setup
set.seed(42)
N <- 5
T_len <- 100

# Independent white noise
X_indep <- matrix(rnorm(N * T_len), nrow = N, ncol = T_len)

# Correlated series
X_corr <- matrix(rnorm(N * T_len), nrow = N, ncol = T_len)
X_corr[2,] <- 0.8 * X_corr[1,] + 0.2 * rnorm(T_len)
X_corr[3,] <- 0.6 * X_corr[1,] + 0.4 * rnorm(T_len)

# List format
X_list <- list(
  s1 = rnorm(T_len),
  s2 = rnorm(T_len),
  s3 = rnorm(T_len)
)

# ============================================================================
# Input validation tests
# ============================================================================

test_that("validate_multivariate_input handles matrix input", {
  result <- validate_multivariate_input(X_indep)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), N)
  expect_equal(ncol(result), T_len)
})

test_that("validate_multivariate_input handles list input", {
  result <- validate_multivariate_input(X_list)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), T_len)
})

test_that("validate_multivariate_input handles data frame", {
  df <- data.frame(x1 = rnorm(100), x2 = rnorm(100), x3 = rnorm(100))
  result <- validate_multivariate_input(df)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 100)
})

test_that("validate_multivariate_input rejects invalid input", {
  expect_error(validate_multivariate_input(matrix(1:10, nrow = 1)), "at least 2")
  expect_error(validate_multivariate_input(matrix(1:4, nrow = 2, ncol = 2)), "at least 3")
  expect_error(validate_multivariate_input(list(rnorm(10))), "at least 2")
  expect_error(validate_multivariate_input(list(rnorm(10), rnorm(15))), "equal length")
})

test_that("validate_multivariate_input handles NA in matrix", {
  # This should NOT error - NA checking is done in main function
  X_na <- X_indep
  X_na[1, 1] <- NA
  result <- validate_multivariate_input(X_na)
  expect_true(anyNA(result))
})

test_that("standardize_multivariate works correctly", {
  X_std <- standardize_multivariate(X_indep)

  # Check each row has mean ~0 and sd ~1
  row_means <- apply(X_std, 1, mean)
  row_sds <- apply(X_std, 1, sd)

  expect_equal(row_means, rep(0, N), tolerance = 1e-10)
  expect_equal(row_sds, rep(1, N), tolerance = 1e-10)
})

test_that("standardize_multivariate handles constant series", {
  X_const <- X_indep
  X_const[1,] <- 5  # Constant series

  expect_warning(X_std <- standardize_multivariate(X_const), "constant")
  expect_equal(X_std[1,], rep(0, T_len))
})

# ============================================================================
# PCA features tests
# ============================================================================

test_that("ts_mv_pca returns correct number of features", {
  result <- ts_mv_pca(X_indep)
  expect_equal(length(result), 15)
  expect_true(all(c("pca_var_pc1", "pca_var_pc2", "pca_effective_rank") %in% names(result)))
})

test_that("ts_mv_pca variance explained sums correctly", {
  result <- ts_mv_pca(X_indep)

  # PC1 + PC2 + PC3 should be less than cumsum_3
  if (!is.na(result$pca_var_pc3)) {
    expect_lte(result$pca_var_pc1 + result$pca_var_pc2 + result$pca_var_pc3,
               result$pca_var_cumsum_3 + 1e-6)
  }
})

test_that("ts_mv_pca detects structure in correlated data", {
  result_indep <- ts_mv_pca(X_indep)
  result_corr <- ts_mv_pca(X_corr)

  # Correlated data should have higher PC1 variance
  expect_gt(result_corr$pca_var_pc1, result_indep$pca_var_pc1)

  # Correlated data should have lower effective rank
  expect_lt(result_corr$pca_effective_rank, result_indep$pca_effective_rank)
})

test_that("ts_mv_pca handles small N correctly", {
  X_small <- X_indep[1:2, ]
  result <- ts_mv_pca(X_small)

  expect_false(is.na(result$pca_var_pc1))
  expect_false(is.na(result$pca_var_pc2))
  expect_true(is.na(result$pca_var_pc3))  # Only 2 series
})

test_that("ts_mv_pca eigenvalue metrics are sensible", {
  result <- ts_mv_pca(X_indep)

  # Effective rank should be between 1 and N
  expect_gte(result$pca_effective_rank, 1)
  expect_lte(result$pca_effective_rank, N)

  # Participation ratio should be between 0 and 1
  expect_gte(result$pca_participation_ratio, 0)
  expect_lte(result$pca_participation_ratio, 1)

  # Entropy should be non-negative
  expect_gte(result$pca_entropy, 0)

  # Stable rank should be between 1 and N
  expect_gte(result$pca_stable_rank, 1)
  expect_lte(result$pca_stable_rank, N)
})

test_that("ts_mv_pca num_pcs thresholds are sensible", {
  result <- ts_mv_pca(X_indep)

  # 95% threshold should need >= 90% threshold
  if (!is.na(result$pca_num_pcs_90pct) && !is.na(result$pca_num_pcs_95pct)) {
    expect_gte(result$pca_num_pcs_95pct, result$pca_num_pcs_90pct)
  }

  # Both should be <= N
  expect_lte(result$pca_num_pcs_90pct, N)
  expect_lte(result$pca_num_pcs_95pct, N)
})

# ============================================================================
# Correlation features tests
# ============================================================================

test_that("ts_mv_correlation returns correct number of features", {
  result <- ts_mv_correlation(X_indep)
  expect_equal(length(result), 15)
  expect_true(all(c("corr_mean", "corr_max", "corr_log_determinant") %in% names(result)))
})

test_that("ts_mv_correlation detects correlations", {
  result_indep <- ts_mv_correlation(X_indep)
  result_corr <- ts_mv_correlation(X_corr)

  # Correlated data should have higher mean correlation
  expect_gt(result_corr$corr_mean, result_indep$corr_mean)

  # Correlated data should have higher max correlation
  expect_gt(result_corr$corr_max, result_indep$corr_max)
})

test_that("ts_mv_correlation bounds are correct", {
  result <- ts_mv_correlation(X_indep)

  # All correlations should be in [-1, 1]
  expect_gte(result$corr_min, -1)
  expect_lte(result$corr_max, 1)

  # Fractions should be in [0, 1]
  expect_gte(result$corr_positive_frac, 0)
  expect_lte(result$corr_positive_frac, 1)
  expect_gte(result$corr_strong_frac, 0)
  expect_lte(result$corr_strong_frac, 1)
  expect_gte(result$corr_weak_frac, 0)
  expect_lte(result$corr_weak_frac, 1)
})

test_that("ts_mv_correlation matrix properties", {
  result <- ts_mv_correlation(X_indep)

  # Spectral radius should be >= 1 (largest eigenvalue of correlation matrix)
  expect_gte(result$corr_spectral_radius, 1)

  # Frobenius norm should be non-negative
  expect_gte(result$corr_frobenius_norm, 0)
})

test_that("ts_mv_correlation handles N=2 correctly", {
  X_two <- X_indep[1:2, ]
  result <- ts_mv_correlation(X_two)

  # Should have 1 correlation value
  expect_false(is.na(result$corr_mean))
  expect_equal(result$corr_mean, result$corr_median)
  expect_equal(result$corr_std, 0)  # Only one value
})

# ============================================================================
# Main function tests
# ============================================================================

test_that("ts_features_multivariate works with matrix input", {
  result <- ts_features_multivariate(X_indep, features = "pca")
  expect_true(is.list(result))
  expect_equal(length(result), 15)
})

test_that("ts_features_multivariate works with list input", {
  result <- ts_features_multivariate(X_list, features = "correlation")
  expect_true(is.list(result))
  expect_equal(length(result), 15)
})

test_that("ts_features_multivariate extracts multiple feature sets", {
  result <- ts_features_multivariate(X_indep, features = c("pca", "correlation"))
  expect_equal(length(result), 15 + 15)  # PCA: 15, Correlation: 15
})

test_that("ts_features_multivariate handles 'all' correctly", {
  result <- ts_features_multivariate(X_indep, features = "all")
  # Phase 1: 15 PCA + 15 correlation
  # Phase 2: 5 covariance + 8 sync + 7 diversity
  # Phase 3: 2 total_correlation
  # Total: 61 features
  expect_equal(length(result), 61)
})

test_that("ts_features_multivariate rejects NA data", {
  X_na <- X_indep
  X_na[1, 1] <- NA
  expect_error(ts_features_multivariate(X_na), "Missing data")
})

test_that("ts_features_multivariate rejects invalid feature sets", {
  expect_error(ts_features_multivariate(X_indep, features = "invalid"), "Invalid feature")
})

test_that("ts_features_multivariate standardization option works", {
  # Without standardization
  X_unscaled <- X_indep * 100  # Scale up
  result_no_std <- ts_features_multivariate(X_unscaled, standardize = FALSE, features = "pca")
  result_std <- ts_features_multivariate(X_unscaled, standardize = TRUE, features = "pca")

  # Results should be identical (PCA on correlation matrix is scale-invariant)
  expect_equal(result_std$pca_var_pc1, result_no_std$pca_var_pc1, tolerance = 0.01)
})

test_that("ts_features_multivariate_df returns data frame", {
  result <- ts_features_multivariate_df(X_indep, features = "pca")
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 15)
})

# ============================================================================
# Integration tests with realistic data
# ============================================================================

test_that("multivariate features work on realistic sensor data", {
  # Simulate 4 vibration sensors with some common mode
  set.seed(123)
  n_sensors <- 4
  n_time <- 200

  # Common vibration signal
  common <- sin(2 * pi * seq(0, 10, length.out = n_time)) + rnorm(n_time, sd = 0.1)

  # Each sensor sees common + individual noise
  sensors <- matrix(0, nrow = n_sensors, ncol = n_time)
  for (i in 1:n_sensors) {
    sensors[i,] <- 0.7 * common + 0.3 * rnorm(n_time)
  }

  # Extract features
  result <- ts_features_multivariate(sensors, features = c("pca", "correlation"))

  # Should detect strong first component (common mode)
  expect_gt(result$pca_var_pc1, 0.5)

  # Should detect high correlations
  expect_gt(result$corr_mean, 0.3)

  # Low effective rank (strong common component)
  expect_lt(result$pca_effective_rank, 3)
})

test_that("multivariate features distinguish independent vs correlated", {
  set.seed(456)

  # Independent series
  X_ind <- matrix(rnorm(3 * 150), nrow = 3)
  feat_ind <- ts_features_multivariate(X_ind, features = "all")

  # Highly correlated series
  X_dep <- matrix(0, nrow = 3, ncol = 150)
  base <- rnorm(150)
  X_dep[1,] <- base + rnorm(150, sd = 0.1)
  X_dep[2,] <- base + rnorm(150, sd = 0.1)
  X_dep[3,] <- base + rnorm(150, sd = 0.1)
  feat_dep <- ts_features_multivariate(X_dep, features = "all")

  # Correlated should have higher PC1 variance
  expect_gt(feat_dep$pca_var_pc1, feat_ind$pca_var_pc1)

  # Correlated should have higher mean correlation
  expect_gt(feat_dep$corr_mean, feat_ind$corr_mean)

  # Correlated should have lower effective rank
  expect_lt(feat_dep$pca_effective_rank, feat_ind$pca_effective_rank)
})
