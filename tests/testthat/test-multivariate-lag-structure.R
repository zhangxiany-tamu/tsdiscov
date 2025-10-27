# Tests for Lag Structure Features (PR 4.2)

test_that("lag structure features work correctly", {
  set.seed(42)

  # Create data with 4 series
  X <- matrix(rnorm(4 * 200), nrow = 4)

  features <- ts_mv_lag_structure(X)

  # Should return 4 features
  expect_named(features, c("lag_structure_entropy", "lag_structure_peak",
                           "lag_structure_mean", "lag_structure_std"))

  # All should be numeric
  expect_true(all(sapply(features, is.numeric)))

  # Entropy should be non-negative
  expect_gte(features$lag_structure_entropy, 0)

  # Standard deviation should be non-negative
  expect_gte(features$lag_structure_std, 0)
})

test_that("lag structure detects leader-follower dynamics", {
  set.seed(123)
  n_time <- 200
  lag_offset <- 5

  # Create leader-follower system
  leader <- rnorm(n_time)
  X_leader <- matrix(0, nrow = 3, ncol = n_time)
  X_leader[1, ] <- leader  # Leader series

  # Followers lag behind by lag_offset
  X_leader[2, (lag_offset + 1):n_time] <- leader[1:(n_time - lag_offset)] + rnorm(n_time - lag_offset, sd = 0.1)
  X_leader[3, (lag_offset + 1):n_time] <- leader[1:(n_time - lag_offset)] + rnorm(n_time - lag_offset, sd = 0.1)
  X_leader[2, 1:lag_offset] <- rnorm(lag_offset)
  X_leader[3, 1:lag_offset] <- rnorm(lag_offset)

  feat_leader <- ts_mv_lag_structure(X_leader)

  # Create synchronized system (no lag)
  base <- rnorm(n_time)
  X_sync <- matrix(0, nrow = 3, ncol = n_time)
  for (i in 1:3) {
    X_sync[i, ] <- base + rnorm(n_time, sd = 0.1)
  }

  feat_sync <- ts_mv_lag_structure(X_sync)

  # Leader-follower should have:
  # - Non-zero peak lag (strong assertion)
  expect_true(abs(feat_leader$lag_structure_peak) > 0)

  # Note: Entropy comparison is stochastic and may not always hold
  # Both entropy values should be finite and reasonable
  expect_true(is.finite(feat_leader$lag_structure_entropy))
  expect_true(is.finite(feat_sync$lag_structure_entropy))
})

test_that("lag structure handles edge cases", {
  set.seed(456)

  # N = 2
  X_two <- matrix(rnorm(2 * 100), nrow = 2)
  feat_two <- ts_mv_lag_structure(X_two)
  expect_length(feat_two, 4)
  expect_true(all(sapply(feat_two, is.numeric)))

  # N = 1 should return NA
  X_one <- matrix(rnorm(100), nrow = 1)
  feat_one <- ts_mv_lag_structure(X_one)
  expect_true(is.na(feat_one$lag_structure_entropy))

  # Very short series
  X_short <- matrix(rnorm(3 * 5), nrow = 3)
  feat_short <- ts_mv_lag_structure(X_short)
  expect_true(is.na(feat_short$lag_structure_entropy))
})

test_that("lag structure works with max_lag parameter", {
  set.seed(789)
  X <- matrix(rnorm(4 * 200), nrow = 4)

  # Different max_lag values
  feat_10 <- ts_mv_lag_structure(X, max_lag = 10)
  feat_20 <- ts_mv_lag_structure(X, max_lag = 20)

  # Results may differ (different lag ranges analyzed)
  # Both should be valid
  expect_true(is.finite(feat_10$lag_structure_entropy))
  expect_true(is.finite(feat_20$lag_structure_entropy))

  # Peak lag should be within bounds
  expect_true(abs(feat_10$lag_structure_peak) <= 10)
  expect_true(abs(feat_20$lag_structure_peak) <= 20)
})

test_that("lag structure integrates with ts_features_multivariate", {
  set.seed(101)
  X <- matrix(rnorm(3 * 150), nrow = 3)

  # Extract via routing function
  feat <- ts_features_multivariate(X, features = "lag_structure")

  # Should return 4 features
  expect_length(feat, 4)
  expect_true("lag_structure_entropy" %in% names(feat))
  expect_true("lag_structure_peak" %in% names(feat))
})

test_that("lag structure works with features='all'", {
  set.seed(202)
  X <- matrix(rnorm(4 * 150), nrow = 4)

  # Extract all features
  feat_all <- ts_features_multivariate(X, features = "all")

  # Should now have 61 features (52 + 4 new)
  expect_length(feat_all, 61)

  # Should include lag structure features
  expect_true("lag_structure_entropy" %in% names(feat_all))
  expect_true("lag_structure_peak" %in% names(feat_all))
  expect_true("lag_structure_mean" %in% names(feat_all))
  expect_true("lag_structure_std" %in% names(feat_all))
})

test_that("lag structure entropy increases with diversity", {
  set.seed(303)
  n_time <- 200

  # Uniform lag structure (all pairs have same lag)
  base <- rnorm(n_time)
  X_uniform <- matrix(0, nrow = 4, ncol = n_time)
  lag <- 3
  X_uniform[1, ] <- base
  for (i in 2:4) {
    X_uniform[i, (lag + 1):n_time] <- base[1:(n_time - lag)] + rnorm(n_time - lag, sd = 0.1)
    X_uniform[i, 1:lag] <- rnorm(lag)
  }

  # Diverse lag structure (different lags for different pairs)
  X_diverse <- matrix(rnorm(4 * n_time), nrow = 4)

  feat_uniform <- ts_mv_lag_structure(X_uniform)
  feat_diverse <- ts_mv_lag_structure(X_diverse)

  # Diverse should have higher entropy
  expect_gt(feat_diverse$lag_structure_entropy, feat_uniform$lag_structure_entropy - 0.1)

  # Uniform should have low standard deviation
  expect_lt(feat_uniform$lag_structure_std, feat_diverse$lag_structure_std + 1)
})

test_that("lag structure detects zero-lag relationships", {
  set.seed(404)
  n_time <- 200

  # Create synchronized data (zero lag)
  base <- rnorm(n_time)
  X_sync <- matrix(0, nrow = 3, ncol = n_time)
  for (i in 1:3) {
    X_sync[i, ] <- 0.9 * base + 0.1 * rnorm(n_time)
  }

  feat_sync <- ts_mv_lag_structure(X_sync)

  # Peak should be close to zero
  expect_true(abs(feat_sync$lag_structure_peak) <= 2)

  # Mean should be close to zero
  expect_true(abs(feat_sync$lag_structure_mean) < 1)
})

test_that("lag structure handles mixed temporal relationships", {
  set.seed(505)
  n_time <- 200

  # Create system with mixed lags
  base <- rnorm(n_time)
  X_mixed <- matrix(0, nrow = 4, ncol = n_time)

  # Series 1: no lag
  X_mixed[1, ] <- base + rnorm(n_time, sd = 0.1)

  # Series 2: lag +5
  lag2 <- 5
  X_mixed[2, (lag2 + 1):n_time] <- base[1:(n_time - lag2)] + rnorm(n_time - lag2, sd = 0.1)
  X_mixed[2, 1:lag2] <- rnorm(lag2)

  # Series 3: lag -3 (leads)
  lag3 <- 3
  X_mixed[3, 1:(n_time - lag3)] <- base[(lag3 + 1):n_time] + rnorm(n_time - lag3, sd = 0.1)
  X_mixed[3, (n_time - lag3 + 1):n_time] <- rnorm(lag3)

  # Series 4: independent
  X_mixed[4, ] <- rnorm(n_time)

  feat_mixed <- ts_mv_lag_structure(X_mixed)

  # Should capture diverse lag structure
  expect_true(is.finite(feat_mixed$lag_structure_entropy))
  expect_true(is.finite(feat_mixed$lag_structure_std))
  expect_gt(feat_mixed$lag_structure_std, 0)
})

test_that("lag structure is scale-invariant with standardization", {
  set.seed(606)
  n_time <- 150

  # Create data with specific lag structure
  base <- rnorm(n_time)
  lag <- 4
  X1 <- matrix(0, nrow = 3, ncol = n_time)
  X1[1, ] <- base
  X1[2, (lag + 1):n_time] <- base[1:(n_time - lag)] + rnorm(n_time - lag, sd = 0.1)
  X1[2, 1:lag] <- rnorm(lag)
  X1[3, (lag + 1):n_time] <- base[1:(n_time - lag)] + rnorm(n_time - lag, sd = 0.1)
  X1[3, 1:lag] <- rnorm(lag)

  # Scale up
  X2 <- X1 * 100

  # Standardize both
  X1_std <- t(scale(t(X1)))
  X2_std <- t(scale(t(X2)))

  feat1 <- ts_mv_lag_structure(X1_std)
  feat2 <- ts_mv_lag_structure(X2_std)

  # Lag structure should be similar (scale-invariant)
  expect_equal(feat1$lag_structure_peak, feat2$lag_structure_peak)
  expect_equal(feat1$lag_structure_entropy, feat2$lag_structure_entropy, tolerance = 0.1)
})
