# Tests for PR 1.2: Parameter Pass-Through
# Verify that parameters are correctly passed through the routing layer

test_that("max_lag parameter affects synchronization features", {
  set.seed(123)
  X <- matrix(rnorm(4 * 100), nrow = 4)

  # Extract with default max_lag = 20
  feat_default <- ts_features_multivariate(X, features = "sync")

  # Extract with custom max_lag = 10
  feat_custom <- ts_features_multivariate(X, features = "sync", max_lag = 10)

  # Results should be different (different lag ranges analyzed)
  expect_false(isTRUE(all.equal(feat_default, feat_custom)))

  # Should still return 8 features
  expect_length(feat_default, 8)
  expect_length(feat_custom, 8)

  # All values should be numeric
  expect_true(all(sapply(feat_default, is.numeric)))
  expect_true(all(sapply(feat_custom, is.numeric)))
})

test_that("max_lag parameter is capped appropriately", {
  set.seed(456)
  # Short series where T/4 < 20
  X_short <- matrix(rnorm(3 * 50), nrow = 3)  # T=50, T/4=12

  # Request max_lag = 100 (will be capped at 12)
  feat <- ts_features_multivariate(X_short, features = "sync", max_lag = 100)

  # Should not error and return 8 features
  expect_length(feat, 8)
  expect_true(all(sapply(feat, is.numeric)))
})

test_that("bins parameter affects synchronization features", {
  set.seed(789)
  X <- matrix(rnorm(4 * 150), nrow = 4)

  # Extract with default bins = 10
  feat_default <- ts_features_multivariate(X, features = "sync")

  # Extract with custom bins = 20
  feat_custom <- ts_features_multivariate(X, features = "sync", bins = 20)

  # Mutual information features should differ (different histogram resolution)
  expect_false(isTRUE(all.equal(feat_default$mutual_info_mean,
                                  feat_custom$mutual_info_mean)))

  # Should still return 8 features
  expect_length(feat_default, 8)
  expect_length(feat_custom, 8)
})

test_that("bins parameter affects diversity features", {
  set.seed(101)
  X <- matrix(rnorm(4 * 150), nrow = 4)

  # Extract with default bins = 10
  feat_default <- ts_features_multivariate(X, features = "diversity")

  # Extract with custom bins = 15
  feat_custom <- ts_features_multivariate(X, features = "diversity", bins = 15)

  # Entropy-based features should differ
  expect_false(isTRUE(all.equal(feat_default$entropy_diversity,
                                  feat_custom$entropy_diversity)))

  # Should still return 7 features
  expect_length(feat_default, 7)
  expect_length(feat_custom, 7)
})

test_that("correlation_method parameter affects correlation features", {
  set.seed(202)
  # Create data with monotonic but nonlinear relationships
  X <- matrix(0, nrow = 3, ncol = 100)
  X[1,] <- rnorm(100)
  X[2,] <- X[1,]^2 + rnorm(100, sd = 0.1)  # Nonlinear relationship
  X[3,] <- rnorm(100)

  # Extract with Pearson (default)
  feat_pearson <- ts_features_multivariate(X, features = "correlation",
                                           correlation_method = "pearson")

  # Extract with Spearman (robust to nonlinear monotonic)
  feat_spearman <- ts_features_multivariate(X, features = "correlation",
                                            correlation_method = "spearman")

  # Results should differ
  expect_false(isTRUE(all.equal(feat_pearson, feat_spearman)))

  # Both should return 15 features
  expect_length(feat_pearson, 15)
  expect_length(feat_spearman, 15)

  # All values should be numeric
  expect_true(all(sapply(feat_pearson, is.numeric)))
  expect_true(all(sapply(feat_spearman, is.numeric)))
})

test_that("kendall correlation method works", {
  set.seed(303)
  X <- matrix(rnorm(3 * 80), nrow = 3)

  # Extract with Kendall
  feat_kendall <- ts_features_multivariate(X, features = "correlation",
                                           correlation_method = "kendall")

  # Should return 15 features
  expect_length(feat_kendall, 15)
  expect_true(all(sapply(feat_kendall, is.numeric)))
})

test_that("invalid correlation method throws error", {
  X <- matrix(rnorm(3 * 100), nrow = 3)

  expect_error(
    ts_features_multivariate(X, features = "correlation",
                             correlation_method = "invalid"),
    "should be one of"
  )
})

test_that("multiple parameters can be passed simultaneously", {
  set.seed(404)
  X <- matrix(rnorm(4 * 150), nrow = 4)

  # Extract with multiple custom parameters
  feat <- ts_features_multivariate(X,
                                   features = c("sync", "correlation", "diversity"),
                                   max_lag = 15,
                                   bins = 12,
                                   correlation_method = "spearman")

  # Should return 15 + 8 + 7 = 30 features
  expect_length(feat, 30)
  expect_true(all(sapply(feat, is.numeric)))
})

test_that("unknown parameters in ... are ignored safely", {
  set.seed(505)
  X <- matrix(rnorm(3 * 100), nrow = 3)

  # Pass unknown parameters - should not error
  expect_silent({
    feat <- ts_features_multivariate(X,
                                     features = "pca",
                                     unknown_param1 = 123,
                                     unknown_param2 = "test")
  })

  # Should still return 15 PCA features
  expect_length(feat, 15)
})

test_that("parameters work with features='all'", {
  set.seed(606)
  X <- matrix(rnorm(4 * 100), nrow = 4)

  # Extract all features with custom parameters
  feat <- ts_features_multivariate(X,
                                   features = "all",
                                   max_lag = 10,
                                   bins = 15,
                                   correlation_method = "spearman")

  # Should return 61 features
  expect_length(expect_length(feat, 61), 61)
  expect_true(all(sapply(feat, is.numeric)))
})

test_that("direct calls to feature functions accept parameters", {
  set.seed(707)
  X <- matrix(rnorm(3 * 100), nrow = 3)

  # Direct call to ts_mv_correlation
  corr_feat <- ts_mv_correlation(X, method = "spearman")
  expect_length(corr_feat, 15)

  # Direct call to ts_mv_synchronization
  sync_feat <- ts_mv_synchronization(X, max_lag = 15, bins = 12)
  expect_length(sync_feat, 8)

  # Direct call to ts_mv_diversity
  div_feat <- ts_mv_diversity(X, bins = 20)
  expect_length(div_feat, 7)
})

test_that("parameter defaults match documented values", {
  set.seed(808)
  X <- matrix(rnorm(3 * 100), nrow = 3)

  # Extract with explicit defaults
  feat_explicit <- ts_features_multivariate(X,
                                            features = "sync",
                                            max_lag = 20,
                                            bins = 10)

  # Extract with implicit defaults
  feat_implicit <- ts_features_multivariate(X, features = "sync")

  # Should be identical
  expect_equal(feat_explicit, feat_implicit)
})
