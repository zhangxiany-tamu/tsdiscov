# Final comprehensive test for 51 multivariate features (after redundancy removal)

test_that("All multivariate features work correctly - 51 total", {
  set.seed(42)
  X <- matrix(rnorm(5 * 200), nrow = 5, ncol = 200)
  
  # Extract all features
  features <- ts_features_multivariate(X, features = "all")
  
  # Should have exactly 51 features
  expect_equal(length(features), 51)
  
  # All should be numeric
  expect_true(all(sapply(features, is.numeric)))
  
  # No duplicated names
  expect_false(any(duplicated(names(features))))
  
  # All names should be non-empty
  expect_true(all(nzchar(names(features))))
})

test_that("Individual feature sets have correct counts", {
  set.seed(42)
  X <- matrix(rnorm(5 * 200), nrow = 5)
  
  # PCA: 15
  pca_feats <- ts_features_multivariate(X, features = "pca")
  expect_equal(length(pca_feats), 15)
  
  # Correlation: 15 (removed corr_trace)
  corr_feats <- ts_features_multivariate(X, features = "correlation")
  expect_equal(length(corr_feats), 15)
  
  # Covariance: 6 (removed cov_nuclear_norm)
  cov_feats <- ts_features_multivariate(X, features = "covariance")
  expect_equal(length(cov_feats), 6)
  
  # Sync: 8 (removed ccf_sync_frac, ccf_zero_lag_mean)
  sync_feats <- ts_features_multivariate(X, features = "sync")
  expect_equal(length(sync_feats), 8)
  
  # Diversity: 7 (removed variance_diversity)
  div_feats <- ts_features_multivariate(X, features = "diversity")
  expect_equal(length(div_feats), 7)
  
  # Total
  expect_equal(15 + 15 + 6 + 8 + 7, 51)
})

test_that("ts_features_all works with all input types", {
  set.seed(42)
  
  # Vector input
  x <- rnorm(200)
  features_vec <- ts_features_all(x)
  expect_true(length(features_vec) > 300)  # Univariate features
  
  # Matrix input
  X <- matrix(rnorm(4 * 200), nrow = 4)
  features_mat <- ts_features_all(X)
  expect_equal(length(features_mat), 51)  # Multivariate features
  
  # List input
  X_list <- list(a = rnorm(150), b = rnorm(150), c = rnorm(150))
  features_list <- ts_features_all(X_list)
  expect_equal(length(features_list), 51)
})

test_that("No redundant features remain", {
  set.seed(42)
  X <- matrix(rnorm(5 * 200), nrow = 5)
  features <- ts_features_multivariate(X, features = "all")
  
  # Removed features should NOT exist
  expect_false("corr_trace" %in% names(features))
  expect_false("cov_nuclear_norm" %in% names(features))
  expect_false("ccf_sync_frac" %in% names(features))
  expect_false("ccf_zero_lag_mean" %in% names(features))
  expect_false("variance_diversity" %in% names(features))
  
  # No spectral features
  expect_false(any(grepl("coherence", names(features))))
  expect_false(any(grepl("spectral_corr", names(features))))
})

test_that("Features work on different data types", {
  set.seed(42)
  
  # Independent signals
  X_indep <- matrix(rnorm(4 * 200), nrow = 4)
  feat_indep <- ts_features_multivariate(X_indep)
  expect_equal(length(feat_indep), 51)
  
  # Correlated signals
  common <- sin(2*pi*seq(0, 10, length.out = 200))
  X_corr <- matrix(0, nrow = 4, ncol = 200)
  for (i in 1:4) X_corr[i,] <- 0.8*common + 0.2*rnorm(200)
  feat_corr <- ts_features_multivariate(X_corr)
  expect_equal(length(feat_corr), 51)
  
  # Should distinguish between independent and correlated
  expect_true(feat_corr$corr_mean > feat_indep$corr_mean)
  expect_true(feat_corr$pca_var_pc1 > feat_indep$pca_var_pc1)
})
