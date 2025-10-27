# Tests for Total Correlation Features (PR 4.1)

test_that("total correlation features work correctly", {
  set.seed(42)

  # Create data with 4 series
  X <- matrix(rnorm(4 * 150), nrow = 4)

  features <- ts_mv_total_correlation(X)

  # Should return 2 features
  expect_named(features, c("total_correlation", "dual_total_correlation"))

  # All should be numeric
  expect_true(all(sapply(features, is.numeric)))

  # Total correlation should be non-negative
  expect_gte(features$total_correlation, 0)

  # Dual total correlation should be non-negative
  expect_gte(features$dual_total_correlation, 0)
})

test_that("total correlation detects independence", {
  set.seed(123)

  # Independent series should have low total correlation
  X_indep <- matrix(rnorm(4 * 200), nrow = 4)
  feat_indep <- ts_mv_total_correlation(X_indep)

  # Create correlated series
  base <- rnorm(200)
  X_corr <- matrix(0, nrow = 4, ncol = 200)
  for (i in 1:4) {
    X_corr[i, ] <- 0.8 * base + 0.2 * rnorm(200)
  }
  feat_corr <- ts_mv_total_correlation(X_corr)

  # Correlated series should have higher total correlation
  expect_gt(feat_corr$total_correlation, feat_indep$total_correlation)

  # Both should be finite
  expect_true(is.finite(feat_indep$total_correlation))
  expect_true(is.finite(feat_corr$total_correlation))
})

test_that("total correlation handles edge cases", {
  set.seed(456)

  # N = 2
  X_two <- matrix(rnorm(2 * 100), nrow = 2)
  feat_two <- ts_mv_total_correlation(X_two)
  expect_length(feat_two, 2)
  expect_true(all(sapply(feat_two, is.numeric)))

  # N = 1 should return NA
  X_one <- matrix(rnorm(100), nrow = 1)
  feat_one <- ts_mv_total_correlation(X_one)
  expect_true(is.na(feat_one$total_correlation))
  expect_true(is.na(feat_one$dual_total_correlation))

  # Very short series
  X_short <- matrix(rnorm(3 * 5), nrow = 3)
  feat_short <- ts_mv_total_correlation(X_short)
  expect_true(is.na(feat_short$total_correlation))
})

test_that("total correlation works with bins parameter", {
  set.seed(789)
  X <- matrix(rnorm(4 * 150), nrow = 4)

  # Different bin sizes
  feat_10 <- ts_mv_total_correlation(X, bins = 10)
  feat_20 <- ts_mv_total_correlation(X, bins = 20)

  # Results should differ (different discretization)
  expect_false(isTRUE(all.equal(feat_10$total_correlation,
                                 feat_20$total_correlation)))

  # Both should be valid
  expect_true(is.finite(feat_10$total_correlation))
  expect_true(is.finite(feat_20$total_correlation))
})

test_that("total correlation integrates with ts_features_multivariate", {
  set.seed(101)
  X <- matrix(rnorm(3 * 150), nrow = 3)

  # Extract via routing function
  feat <- ts_features_multivariate(X, features = "total_correlation")

  # Should return 2 features
  expect_length(feat, 2)
  expect_true("total_correlation" %in% names(feat))
  expect_true("dual_total_correlation" %in% names(feat))
})

test_that("total correlation works with features='all'", {
  set.seed(202)
  X <- matrix(rnorm(4 * 150), nrow = 4)

  # Extract all features
  feat_all <- ts_features_multivariate(X, features = "all")

  # Should now have 61 features (50 + 2 new)
  expect_length(feat_all, 61)

  # Should include total correlation features
  expect_true("total_correlation" %in% names(feat_all))
  expect_true("dual_total_correlation" %in% names(feat_all))
})

test_that("total correlation handles constant series", {
  set.seed(303)

  # Create data with one constant series
  X <- matrix(rnorm(3 * 100), nrow = 3)
  X[1, ] <- 5  # Constant

  feat <- ts_mv_total_correlation(X)

  # Should handle gracefully (may return NA or valid value)
  expect_true(is.numeric(feat$total_correlation))
  expect_true(is.numeric(feat$dual_total_correlation))
})

test_that("joint entropy computation works", {
  set.seed(404)

  # Test internal function
  X <- matrix(rnorm(3 * 100), nrow = 3)

  # Compute joint entropy
  joint_ent <- compute_joint_entropy(X, bins = 10)

  # Should be positive
  expect_gt(joint_ent, 0)
  expect_true(is.finite(joint_ent))

  # Single series should match marginal entropy
  X_one <- matrix(rnorm(100), nrow = 1)
  joint_one <- compute_joint_entropy(X_one, bins = 10)
  marginal_one <- cpp_mv_shannon_entropy(X_one[1, ], bins = 10)

  expect_equal(joint_one, marginal_one, tolerance = 0.01)
})

test_that("total correlation increases with dependency strength", {
  set.seed(505)
  n_time <- 200

  # Low dependency
  base <- rnorm(n_time)
  X_low <- matrix(0, nrow = 3, ncol = n_time)
  for (i in 1:3) {
    X_low[i, ] <- 0.3 * base + 0.7 * rnorm(n_time)
  }

  # High dependency
  X_high <- matrix(0, nrow = 3, ncol = n_time)
  for (i in 1:3) {
    X_high[i, ] <- 0.9 * base + 0.1 * rnorm(n_time)
  }

  feat_low <- ts_mv_total_correlation(X_low)
  feat_high <- ts_mv_total_correlation(X_high)

  # Higher dependency should give higher total correlation
  expect_gt(feat_high$total_correlation, feat_low$total_correlation)
})

test_that("total correlation handles large N efficiently", {
  set.seed(606)

  # 10 series (should use sampling internally)
  X_large <- matrix(rnorm(10 * 300), nrow = 10)

  # Should complete without error
  feat <- ts_mv_total_correlation(X_large, bins = 10)

  expect_length(feat, 2)
  expect_true(is.finite(feat$total_correlation))
  expect_true(is.finite(feat$dual_total_correlation))
})
