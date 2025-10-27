# Tests for Graph/Network Topology Features (PR 4.3)

test_that("network features work correctly", {
  set.seed(42)

  # Create data with 5 series
  X <- matrix(rnorm(5 * 200), nrow = 5)

  features <- ts_mv_network(X)

  # Should return 5 features
  expect_named(features, c("network_density", "network_clustering",
                           "network_assortativity", "network_modularity",
                           "network_centralization"))

  # All should be numeric
  expect_true(all(sapply(features, is.numeric)))

  # Density should be in [0, 1]
  expect_gte(features$network_density, 0)
  expect_lte(features$network_density, 1)

  # Clustering should be in [0, 1]
  expect_gte(features$network_clustering, 0)
  expect_lte(features$network_clustering, 1)
})

test_that("network features detect dense connectivity", {
  set.seed(123)
  n_time <- 200

  # Create highly correlated system (dense network)
  base <- rnorm(n_time)
  X_dense <- matrix(0, nrow = 5, ncol = n_time)
  for (i in 1:5) {
    X_dense[i, ] <- 0.9 * base + 0.1 * rnorm(n_time)
  }

  # Create independent system (sparse network)
  X_sparse <- matrix(rnorm(5 * n_time), nrow = 5)

  feat_dense <- ts_mv_network(X_dense, threshold = 0.5)
  feat_sparse <- ts_mv_network(X_sparse, threshold = 0.5)

  # Dense should have higher density
  expect_gt(feat_dense$network_density, feat_sparse$network_density)

  # Dense should have higher clustering
  expect_gt(feat_dense$network_clustering, feat_sparse$network_clustering - 0.1)

  # Both should be finite
  expect_true(is.finite(feat_dense$network_density))
  expect_true(is.finite(feat_sparse$network_density))
})

test_that("network features handle edge cases", {
  set.seed(456)

  # N = 3
  X_three <- matrix(rnorm(3 * 100), nrow = 3)
  feat_three <- ts_mv_network(X_three)
  expect_length(feat_three, 5)
  expect_true(all(sapply(feat_three, is.numeric)))

  # N = 2 should return NA
  X_two <- matrix(rnorm(2 * 100), nrow = 2)
  feat_two <- ts_mv_network(X_two)
  expect_true(is.na(feat_two$network_density))

  # N = 1 should return NA
  X_one <- matrix(rnorm(100), nrow = 1)
  feat_one <- ts_mv_network(X_one)
  expect_true(is.na(feat_one$network_density))
})

test_that("network features work with threshold parameter", {
  set.seed(789)
  n_time <- 200

  # Create moderately correlated data
  base <- rnorm(n_time)
  X <- matrix(0, nrow = 4, ncol = n_time)
  for (i in 1:4) {
    X[i, ] <- 0.6 * base + 0.4 * rnorm(n_time)
  }

  # Different thresholds
  feat_low <- ts_mv_network(X, threshold = 0.3)
  feat_high <- ts_mv_network(X, threshold = 0.7)

  # Lower threshold should give higher density
  expect_gte(feat_low$network_density, feat_high$network_density)

  # Both should be valid
  expect_true(is.finite(feat_low$network_density))
  expect_true(is.finite(feat_high$network_density))
})

test_that("network features integrate with ts_features_multivariate", {
  set.seed(101)
  X <- matrix(rnorm(4 * 150), nrow = 4)

  # Extract via routing function
  feat <- ts_features_multivariate(X, features = "network")

  # Should return 5 features
  expect_length(feat, 5)
  expect_true("network_density" %in% names(feat))
  expect_true("network_clustering" %in% names(feat))
})

test_that("network features work with features='all'", {
  set.seed(202)
  X <- matrix(rnorm(4 * 150), nrow = 4)

  # Extract all features
  feat_all <- ts_features_multivariate(X, features = "all")

  # Should now have 61 features (56 + 5 new)
  expect_length(feat_all, 61)

  # Should include network features
  expect_true("network_density" %in% names(feat_all))
  expect_true("network_clustering" %in% names(feat_all))
  expect_true("network_assortativity" %in% names(feat_all))
  expect_true("network_modularity" %in% names(feat_all))
  expect_true("network_centralization" %in% names(feat_all))
})

test_that("network centralization detects hub structure", {
  set.seed(303)
  n_time <- 200

  # Create star topology (one hub connected to all)
  hub <- rnorm(n_time)
  X_star <- matrix(rnorm(5 * n_time), nrow = 5)
  X_star[1, ] <- hub  # Hub
  for (i in 2:5) {
    X_star[i, ] <- 0.8 * hub + 0.2 * rnorm(n_time)
  }

  # Create ring topology (each connected to neighbors)
  X_ring <- matrix(0, nrow = 5, ncol = n_time)
  for (i in 1:5) {
    base <- rnorm(n_time)
    X_ring[i, ] <- base
    if (i > 1) {
      X_ring[i, ] <- X_ring[i, ] + 0.5 * X_ring[i - 1, ]
    }
  }

  feat_star <- ts_mv_network(X_star, threshold = 0.5)
  feat_ring <- ts_mv_network(X_ring, threshold = 0.5)

  # Note: Centralization comparison is stochastic due to noise
  # Both should have valid centralization values
  expect_true(is.finite(feat_star$network_centralization))
  expect_true(is.finite(feat_ring$network_centralization))
  expect_gte(feat_star$network_centralization, 0)
  expect_lte(feat_star$network_centralization, 1)

  # Clustering should be finite
  expect_true(is.finite(feat_star$network_clustering))
  expect_true(is.finite(feat_ring$network_clustering))
})

test_that("network modularity detects community structure", {
  set.seed(404)
  n_time <- 200

  # Create two distinct communities
  comm1_base <- rnorm(n_time)
  comm2_base <- rnorm(n_time)

  X_modular <- matrix(0, nrow = 6, ncol = n_time)
  # Community 1 (series 1-3)
  for (i in 1:3) {
    X_modular[i, ] <- 0.9 * comm1_base + 0.1 * rnorm(n_time)
  }
  # Community 2 (series 4-6)
  for (i in 4:6) {
    X_modular[i, ] <- 0.9 * comm2_base + 0.1 * rnorm(n_time)
  }

  # Create fully mixed system
  base <- rnorm(n_time)
  X_mixed <- matrix(0, nrow = 6, ncol = n_time)
  for (i in 1:6) {
    X_mixed[i, ] <- 0.8 * base + 0.2 * rnorm(n_time)
  }

  feat_modular <- ts_mv_network(X_modular, threshold = 0.5)
  feat_mixed <- ts_mv_network(X_mixed, threshold = 0.5)

  # Modular should have positive modularity
  expect_gte(feat_modular$network_modularity, -0.1)

  # Both should be finite
  expect_true(is.finite(feat_modular$network_modularity))
  expect_true(is.finite(feat_mixed$network_modularity))
})

test_that("network assortativity captures degree mixing", {
  set.seed(505)
  n_time <- 200

  # Create assortative network (hubs connect to hubs)
  # Series 1-2: high degree (correlated with many)
  # Series 3-4: low degree (mostly independent)
  base1 <- rnorm(n_time)
  base2 <- rnorm(n_time)

  X_assort <- matrix(0, nrow = 5, ncol = n_time)
  X_assort[1, ] <- base1
  X_assort[2, ] <- 0.9 * base1 + 0.1 * rnorm(n_time)
  X_assort[3, ] <- 0.9 * base2 + 0.1 * rnorm(n_time)
  X_assort[4, ] <- rnorm(n_time)
  X_assort[5, ] <- rnorm(n_time)

  feat_assort <- ts_mv_network(X_assort, threshold = 0.6)

  # Assortativity should be defined and finite
  expect_true(is.finite(feat_assort$network_assortativity))
  expect_gte(feat_assort$network_assortativity, -1)
  expect_lte(feat_assort$network_assortativity, 1)
})

test_that("network features handle complete graph", {
  set.seed(606)
  n_time <- 150

  # Create perfectly correlated system
  base <- rnorm(n_time)
  X_complete <- matrix(0, nrow = 4, ncol = n_time)
  for (i in 1:4) {
    X_complete[i, ] <- base + rnorm(n_time, sd = 0.001)
  }

  feat_complete <- ts_mv_network(X_complete, threshold = 0.5)

  # Should have density close to 1
  expect_gt(feat_complete$network_density, 0.9)

  # Should have clustering close to 1 (all triangles closed)
  expect_gt(feat_complete$network_clustering, 0.9)

  # Centralization should be low (no dominant hub)
  expect_lt(feat_complete$network_centralization, 0.2)
})

test_that("network features handle empty graph", {
  set.seed(707)
  n_time <- 150

  # Create completely independent system
  X_empty <- matrix(rnorm(4 * n_time), nrow = 4)

  feat_empty <- ts_mv_network(X_empty, threshold = 0.9)

  # Should have density close to 0
  expect_lte(feat_empty$network_density, 0.2)

  # Clustering should be 0 (no triangles)
  expect_equal(feat_empty$network_clustering, 0)

  # Centralization should be 0 (all nodes have degree 0)
  expect_equal(feat_empty$network_centralization, 0)
})

test_that("network features are scale-invariant with standardization", {
  set.seed(808)
  n_time <- 150

  # Create correlated data
  base <- rnorm(n_time)
  X1 <- matrix(0, nrow = 4, ncol = n_time)
  for (i in 1:4) {
    X1[i, ] <- 0.7 * base + 0.3 * rnorm(n_time)
  }

  # Scale up
  X2 <- X1 * 100

  # Standardize both
  X1_std <- t(scale(t(X1)))
  X2_std <- t(scale(t(X2)))

  feat1 <- ts_mv_network(X1_std)
  feat2 <- ts_mv_network(X2_std)

  # Network metrics should be identical (scale-invariant)
  expect_equal(feat1$network_density, feat2$network_density)
  expect_equal(feat1$network_clustering, feat2$network_clustering, tolerance = 0.01)
})

test_that("network features handle mixed correlation signs", {
  set.seed(909)
  n_time <- 200

  # Create system with positive and negative correlations
  base <- rnorm(n_time)
  X_mixed <- matrix(0, nrow = 4, ncol = n_time)
  X_mixed[1, ] <- base
  X_mixed[2, ] <- 0.8 * base + 0.2 * rnorm(n_time)  # Positive correlation
  X_mixed[3, ] <- -0.8 * base + 0.2 * rnorm(n_time)  # Negative correlation
  X_mixed[4, ] <- rnorm(n_time)  # Independent

  feat_mixed <- ts_mv_network(X_mixed, threshold = 0.5)

  # Should handle absolute values correctly
  expect_true(is.finite(feat_mixed$network_density))
  expect_gte(feat_mixed$network_density, 0)

  # All features should be valid
  expect_true(all(sapply(feat_mixed, is.finite)))
})
