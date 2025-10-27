test_that("ts_features_all works with univariate input (vector)", {
  set.seed(42)
  x <- rnorm(200)

  # Auto detection
  features <- ts_features_all(x)
  expect_true(is.list(features))
  expect_true(length(features) > 300)  # Should be ~352 univariate features

  # Explicit univariate
  features2 <- ts_features_all(x, feature_type = "univariate")
  expect_equal(length(features), length(features2))
  expect_equal(names(features), names(features2))

  # Both (should default to univariate for single series)
  expect_warning(
    features3 <- ts_features_all(x, feature_type = "both"),
    "Single series"
  )
  expect_equal(length(features), length(features3))
})

test_that("ts_features_all works with multivariate input (matrix)", {
  set.seed(42)
  X <- matrix(rnorm(5 * 200), nrow = 5, ncol = 200)

  # Auto detection (should return multivariate)
  features <- ts_features_all(X)
  expect_true(is.list(features))
  expect_equal(length(features), 61)  # 61 multivariate features

  # Explicit multivariate
  features2 <- ts_features_all(X, feature_type = "multivariate")
  expect_equal(features, features2)
})

test_that("ts_features_all works with list input", {
  set.seed(42)
  X <- list(
    series1 = rnorm(150),
    series2 = rnorm(150),
    series3 = rnorm(150)
  )

  # Auto detection
  features <- ts_features_all(X)
  expect_true(is.list(features))
  expect_equal(length(features), 61)
})

test_that("ts_features_all feature_type selection works", {
  set.seed(42)
  X <- matrix(rnorm(4 * 200), nrow = 4)

  # Univariate only with per_series
  features_univ <- ts_features_all(X,
                                    feature_type = "univariate",
                                    univariate_summary = "per_series")
  expect_true(all(grepl("^series\\d+_", names(features_univ))))
  expect_true(length(features_univ) > 1000)  # 4 series × ~352 features

  # Multivariate only
  features_mult <- ts_features_all(X, feature_type = "multivariate")
  expect_equal(length(features_mult), 61)
  expect_true(all(!grepl("^series\\d+_", names(features_mult))))

  # Both
  features_both <- ts_features_all(X,
                                    feature_type = "both",
                                    univariate_summary = "aggregate")
  expect_true(length(features_both) > 50)  # Multivariate + aggregated univariate
  expect_true(any(grepl("_mean$", names(features_both))))  # Aggregated features
  expect_true("pca_var_pc1" %in% names(features_both))     # Multivariate features
})

test_that("ts_features_all univariate_summary options work", {
  set.seed(42)
  X <- matrix(rnorm(3 * 150), nrow = 3)

  # None
  features_none <- ts_features_all(X,
                                    feature_type = "both",
                                    univariate_summary = "none")
  expect_equal(length(features_none), 61)  # Only multivariate

  # Aggregate
  features_agg <- ts_features_all(X,
                                   feature_type = "both",
                                   univariate_summary = "aggregate")
  expect_true(length(features_agg) > 50)
  expect_true(any(grepl("_mean$", names(features_agg))))
  expect_true(any(grepl("_max$", names(features_agg))))
  expect_true(any(grepl("_min$", names(features_agg))))
  expect_true(any(grepl("_sd$", names(features_agg))))

  # Per series
  features_per <- ts_features_all(X,
                                   feature_type = "both",
                                   univariate_summary = "per_series")
  expect_true(length(features_per) > 1000)  # 50 + 3×352
  expect_true(any(grepl("^series1_", names(features_per))))
  expect_true(any(grepl("^series2_", names(features_per))))
  expect_true(any(grepl("^series3_", names(features_per))))
})

test_that("ts_features_all multivariate_sets selection works", {
  set.seed(42)
  X <- matrix(rnorm(5 * 200), nrow = 5)

  # Specific sets
  features_sync <- ts_features_all(X,
                                    feature_type = "multivariate",
                                    multivariate_sets = "sync")
  expect_equal(length(features_sync), 8)

  features_pca <- ts_features_all(X,
                                   feature_type = "multivariate",
                                   multivariate_sets = "pca")
  expect_equal(length(features_pca), 15)

  # Multiple sets
  features_combined <- ts_features_all(X,
                                        feature_type = "multivariate",
                                        multivariate_sets = c("sync", "diversity"))
  expect_equal(length(features_combined), 15)  # 8 + 7
})

test_that("ts_features_all error handling works", {
  set.seed(42)
  x <- rnorm(200)

  # Cannot get multivariate from single series
  expect_error(
    ts_features_all(x, feature_type = "multivariate"),
    "Cannot extract multivariate features from single time series"
  )
})

test_that("ts_features_all_df returns data frame", {
  set.seed(42)
  X <- matrix(rnorm(4 * 150), nrow = 4)

  df <- ts_features_all_df(X)
  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 1)
  expect_equal(ncol(df), 61)

  # With both features
  df_both <- ts_features_all_df(X,
                                 feature_type = "both",
                                 univariate_summary = "aggregate")
  expect_true(is.data.frame(df_both))
  expect_equal(nrow(df_both), 1)
  expect_true(ncol(df_both) > 50)
})

test_that("aggregate_univariate_features works correctly", {
  set.seed(42)
  X <- matrix(rnorm(3 * 200), nrow = 3)

  agg_features <- aggregate_univariate_features(X)

  # Should have mean, max, min, sd for each univariate feature
  expect_true(length(agg_features) > 300 * 4)  # At least 4 stats per feature

  # Check specific aggregations exist
  expect_true(any(grepl("_mean$", names(agg_features))))
  expect_true(any(grepl("_max$", names(agg_features))))
  expect_true(any(grepl("_min$", names(agg_features))))
  expect_true(any(grepl("_sd$", names(agg_features))))

  # All values should be numeric
  expect_true(all(sapply(agg_features, is.numeric)))
})

test_that("ts_features_all handles edge cases", {
  set.seed(42)

  # Single series as 1-row matrix
  X_single <- matrix(rnorm(200), nrow = 1)
  features <- ts_features_all(X_single)
  expect_true(length(features) > 300)  # Should treat as univariate

  # Two series (minimum for multivariate)
  X_two <- matrix(rnorm(2 * 150), nrow = 2)
  features_mv <- ts_features_all(X_two, feature_type = "multivariate")
  expect_equal(length(features_mv), 61)
})

test_that("ts_features_all produces consistent feature names", {
  set.seed(42)
  X <- matrix(rnorm(4 * 200), nrow = 4)

  # Multivariate feature names
  features_mv <- ts_features_all(X, feature_type = "multivariate")
  expect_true(all(nzchar(names(features_mv))))
  expect_false(any(duplicated(names(features_mv))))

  # Aggregated feature names
  features_agg <- ts_features_all(X,
                                   feature_type = "both",
                                   univariate_summary = "aggregate")
  expect_true(all(nzchar(names(features_agg))))
  expect_false(any(duplicated(names(features_agg))))

  # Per-series feature names
  features_per <- ts_features_all(X,
                                   feature_type = "both",
                                   univariate_summary = "per_series")
  expect_true(all(nzchar(names(features_per))))
  expect_false(any(duplicated(names(features_per))))
})
