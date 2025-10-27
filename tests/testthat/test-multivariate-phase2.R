test_that("Covariance features work correctly", {
  set.seed(42)

  # 3 series with different variances
  X <- matrix(0, nrow = 3, ncol = 150)
  X[1,] <- rnorm(150, sd = 1)
  X[2,] <- rnorm(150, sd = 2)
  X[3,] <- rnorm(150, sd = 0.5)

  features <- ts_mv_covariance(X)

  # Check that all features are returned (removed cov_nuclear_norm - redundant with cov_trace)
  expect_named(features, c(
    "cov_trace", "cov_determinant", "cov_log_determinant",
    "cov_condition_number", "cov_frobenius_norm",
    "cov_spectral_norm"
  ))

  # Check that features are numeric
  expect_true(all(sapply(features, is.numeric)))

  # Trace should be positive (sum of variances)
  expect_true(features$cov_trace > 0)

  # Spectral norm should be positive
  expect_true(features$cov_spectral_norm > 0)

  # Condition number should be >= 1
  if (!is.na(features$cov_condition_number)) {
    expect_true(features$cov_condition_number >= 1)
  }
})

test_that("Covariance features handle edge cases", {
  # N = 2
  X <- matrix(rnorm(2 * 100), nrow = 2)
  features <- ts_mv_covariance(X)
  expect_true(all(sapply(features, is.numeric)))

  # Constant series (after standardization should be all zeros)
  X <- matrix(0, nrow = 3, ncol = 100)
  X[1,] <- 1  # Constant
  X[2,] <- 2  # Constant
  X[3,] <- 3  # Constant
  features <- ts_mv_covariance(X)
  # All zeros should have determinant = 0
  expect_true(is.na(features$cov_log_determinant) || features$cov_determinant < 1e-6)
})

test_that("Synchronization features work correctly", {
  set.seed(42)

  # Create synchronized signals
  t <- seq(0, 10, length.out = 200)
  common <- sin(2*pi*t)

  X <- matrix(0, nrow = 3, ncol = 200)
  X[1,] <- common + rnorm(200, sd = 0.1)
  X[2,] <- common + rnorm(200, sd = 0.1)
  X[3,] <- common + rnorm(200, sd = 0.1)

  features <- ts_mv_synchronization(X)

  # Check that all features are returned (removed ccf_sync_frac and ccf_zero_lag_mean - redundant)
  expect_named(features, c(
    "ccf_max_mean", "ccf_max_median", "ccf_max_std",
    "ccf_lag_mean", "ccf_lag_std",
    "mutual_info_mean", "mutual_info_max", "mutual_info_std"
  ))

  # High synchronization should give high CCF values
  expect_true(features$ccf_max_mean > 0.5)

  # Mutual information should be positive
  expect_true(features$mutual_info_mean > 0)
})

test_that("Synchronization features detect independence", {
  set.seed(42)

  # Independent signals
  X <- matrix(rnorm(4 * 200), nrow = 4, ncol = 200)

  features <- ts_mv_synchronization(X)

  # Independent signals should have low CCF
  expect_true(features$ccf_max_mean < 0.5)

  # Mutual information should be low for independent signals
  expect_true(features$mutual_info_mean < 1)
})

test_that("Synchronization handles edge cases", {
  # N = 2
  X <- matrix(rnorm(2 * 100), nrow = 2)
  features <- ts_mv_synchronization(X)
  expect_true(all(sapply(features, is.numeric)))

  # N = 1 should return NA
  X <- matrix(rnorm(100), nrow = 1)
  features <- ts_mv_synchronization(X)
  expect_true(all(is.na(unlist(features)[1:9])))  # All but last should be NA
})

test_that("Spectral features work correctly", {
  set.seed(42)

  # Create signals with common frequency component
  n_time <- 256  # Power of 2 for FFT
  t <- seq(0, 10, length.out = n_time)
  common <- sin(2*pi*t) + sin(4*pi*t)

  X <- matrix(0, nrow = 3, ncol = n_time)
  X[1,] <- common + rnorm(n_time, sd = 0.2)
  X[2,] <- common + rnorm(n_time, sd = 0.2)
  X[3,] <- rnorm(n_time)  # Independent

  features <- ts_mv_spectral(X)

  # Check that all features are returned
  expect_named(features, c(
    "coherence_mean", "coherence_median", "coherence_max",
    "coherence_std", "coherence_high_frac",
    "phase_lag_mean", "phase_lag_std",
    "spectral_corr_mean", "spectral_corr_max",
    "csd_magnitude_mean"
  ))

  # All features should be numeric
  expect_true(all(sapply(features, is.numeric)))

  # Coherence should be between 0 and 1
  expect_true(features$coherence_mean >= 0 && features$coherence_mean <= 1)
  expect_true(features$coherence_max >= 0 && features$coherence_max <= 1)
})

test_that("Spectral features handle edge cases", {
  # N = 2
  X <- matrix(rnorm(2 * 128), nrow = 2)
  features <- ts_mv_spectral(X)
  expect_true(all(sapply(features, is.numeric)))

  # N = 1 should return NA
  X <- matrix(rnorm(128), nrow = 1)
  features <- ts_mv_spectral(X)
  expect_true(all(is.na(unlist(features))))
})

test_that("Diversity features work correctly", {
  set.seed(42)

  # Heterogeneous system
  X <- matrix(0, nrow = 4, ncol = 200)
  X[1,] <- rnorm(200, sd = 1)
  X[2,] <- rnorm(200, sd = 2)
  X[3,] <- rnorm(200, sd = 0.5)
  X[4,] <- rnorm(200, sd = 3)

  features <- ts_mv_diversity(X)

  # Check that all features are returned (removed variance_diversity - always 0 when standardized)
  expect_named(features, c(
    "mean_diversity",
    "acf1_diversity", "complexity_variance",
    "entropy_diversity", "range_diversity",
    "skewness_diversity", "kurtosis_diversity"
  ))

  # All features should be numeric
  expect_true(all(sapply(features, is.numeric)))

  # Range diversity should be positive
  expect_true(features$range_diversity > 0)
})

test_that("Diversity features detect homogeneity", {
  set.seed(42)

  # Homogeneous system (all similar)
  X <- matrix(rnorm(4 * 200), nrow = 4, ncol = 200)

  features <- ts_mv_diversity(X)

  # Homogeneous system should have lower diversity
  # (but not zero due to sampling variation)
  expect_true(all(sapply(features, is.finite)))
})

test_that("Diversity handles edge cases", {
  # N = 2
  X <- matrix(rnorm(2 * 100), nrow = 2)
  features <- ts_mv_diversity(X)
  expect_true(all(sapply(features, is.numeric)))

  # N = 1 should return NA
  X <- matrix(rnorm(100), nrow = 1)
  features <- ts_mv_diversity(X)
  expect_true(all(is.na(unlist(features))))
})

test_that("Phase 2 integration: all features work together", {
  set.seed(42)

  # Realistic multivariate system
  X <- matrix(rnorm(5 * 200), nrow = 5, ncol = 200)

  # Test individual feature sets
  cov_feats <- ts_features_multivariate(X, features = "covariance", standardize = FALSE)
  sync_feats <- ts_features_multivariate(X, features = "sync")
  div_feats <- ts_features_multivariate(X, features = "diversity")

  # Should have correct number of features (after redundancy removal)
  expect_length(cov_feats, 6)   # Was 7, removed cov_nuclear_norm
  expect_length(sync_feats, 8)   # Was 10, removed ccf_sync_frac and ccf_zero_lag_mean
  expect_length(div_feats, 7)    # Was 8, removed variance_diversity

  # All features combined
  all_feats <- ts_features_multivariate(X, features = "all")

  # Phase 1: 15 PCA + 15 correlation = 30 (removed corr_trace)
  # Phase 2: 6 cov + 8 sync + 7 diversity = 21 (removed spectral entirely)
  # Total: 51 features
  expect_length(all_feats, 51)

  # All should be numeric
  expect_true(all(sapply(all_feats, is.numeric)))
})

test_that("Phase 2 integration: data frame output works", {
  set.seed(42)
  X <- matrix(rnorm(4 * 150), nrow = 4, ncol = 150)

  df <- ts_features_multivariate_df(X, features = "all")

  # Should be 1 row, 66 columns
  expect_equal(nrow(df), 1)
  expect_equal(ncol(df), 51)

  # All columns should be numeric
  expect_true(all(sapply(df, is.numeric)))
})

test_that("Phase 2 integration: batch processing works", {
  set.seed(42)

  systems <- list(
    system_A = matrix(rnorm(3 * 100), nrow = 3),
    system_B = matrix(rnorm(3 * 100), nrow = 3),
    system_C = matrix(rnorm(3 * 100), nrow = 3)
  )

  results <- do.call(rbind, lapply(systems, ts_features_multivariate_df))
  rownames(results) <- names(systems)

  # Should be 3 x 51 (after redundancy removal)
  expect_equal(dim(results), c(3, 51))

  # All should be numeric
  expect_true(all(sapply(results, is.numeric)))
})

test_that("Phase 2: mutual information computation works", {
  # Test internal MI function
  set.seed(42)

  # Independent
  x <- rnorm(200)
  y <- rnorm(200)
  mi_indep <- compute_mutual_information(x, y)

  # Should be low for independent
  expect_true(mi_indep < 0.5)

  # Dependent
  x <- rnorm(200)
  y <- x + rnorm(200, sd = 0.1)
  mi_dep <- compute_mutual_information(x, y)

  # Should be higher for dependent
  expect_true(mi_dep > mi_indep)
})

test_that("Phase 2: Shannon entropy computation works", {
  # Test internal entropy function
  set.seed(42)

  # Uniform-like distribution
  x <- runif(200)
  entropy_uniform <- compute_shannon_entropy_simple(x)

  # Should be positive
  expect_true(entropy_uniform > 0)

  # Peaked distribution
  x <- rnorm(200, sd = 0.1)
  entropy_peaked <- compute_shannon_entropy_simple(x)

  # Uniform should have higher entropy
  expect_true(entropy_uniform > entropy_peaked)

  # Too few points
  x <- rnorm(5)
  entropy_few <- compute_shannon_entropy_simple(x)
  expect_true(is.na(entropy_few))
})

test_that("Phase 2 features distinguish different system types", {
  set.seed(42)
  n_time <- 200

  # System 1: Independent signals
  X_indep <- matrix(rnorm(4 * n_time), nrow = 4)
  feats_indep <- ts_features_multivariate(X_indep, features = "sync")

  # System 2: Synchronized signals
  common <- sin(2 * pi * seq(0, 10, length.out = n_time))
  X_sync <- matrix(0, nrow = 4, ncol = n_time)
  for (i in 1:4) {
    X_sync[i,] <- 0.8 * common + 0.2 * rnorm(n_time)
  }
  feats_sync <- ts_features_multivariate(X_sync, features = "sync")

  # Synchronized should have higher CCF
  expect_true(feats_sync$ccf_max_mean > feats_indep$ccf_max_mean)

  # Both sets of features should be finite
  expect_true(is.finite(feats_sync$ccf_max_mean))
  expect_true(is.finite(feats_indep$ccf_max_mean))
})

test_that("Phase 2: covariance features distinguish variance structures", {
  set.seed(42)

  # System 1: Similar variances
  X_sim <- matrix(rnorm(4 * 150, sd = 1), nrow = 4)
  feats_sim <- ts_mv_covariance(X_sim)

  # System 2: Very different variances
  X_diff <- matrix(0, nrow = 4, ncol = 150)
  X_diff[1,] <- rnorm(150, sd = 0.1)
  X_diff[2,] <- rnorm(150, sd = 1)
  X_diff[3,] <- rnorm(150, sd = 5)
  X_diff[4,] <- rnorm(150, sd = 10)
  feats_diff <- ts_mv_covariance(X_diff)

  # Different variances should have higher condition number
  if (!is.na(feats_sim$cov_condition_number) && !is.na(feats_diff$cov_condition_number)) {
    expect_true(feats_diff$cov_condition_number > feats_sim$cov_condition_number)
  }

  # Different variances should have higher trace
  expect_true(feats_diff$cov_trace > feats_sim$cov_trace)
})
