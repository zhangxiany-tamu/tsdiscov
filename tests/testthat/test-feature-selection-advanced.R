# Advanced Feature Selection Tests
# Tests for FWER methods and different target types

library(testthat)
library(tsdiscov)

# Helper function to generate test data
generate_test_data <- function(n_per_class = 15, n_obs = 100) {
  set.seed(123)
  # Class 0: White noise
  class0 <- lapply(1:n_per_class, function(i) rnorm(n_obs))
  # Class 1: Random walk
  class1 <- lapply(1:n_per_class, function(i) cumsum(rnorm(n_obs)))

  # Suppress HoltWinters warnings (expected for white noise)
  X <- suppressWarnings({
    rbind(
      do.call(rbind, lapply(class0, ts_features)),
      do.call(rbind, lapply(class1, ts_features))
    )
  })
  y <- c(rep(0, n_per_class), rep(1, n_per_class))

  list(X = X, y = y)
}

# ==============================================================================
# Test FWER Control Methods
# ==============================================================================

test_that("FWER: Bonferroni correction works (most conservative)", {
  data <- generate_test_data()

  # Bonferroni is most conservative (robust to any dependence)
  # Suppress warnings about NA features (expected)
  result <- suppressWarnings({
    select_features(
      data$X, data$y,
      fdr_level = 0.05,
      correction_method = "bonferroni"
    )
  })

  expect_s3_class(result, "tsdiscov_feature_selection")
  expect_gte(result$n_features_selected, 1)
  expect_lte(result$n_features_selected, result$n_features_original)

  # Check correction method is stored
  expect_true("correction_method" %in% names(result$relevance_table))
  expect_equal(unique(result$relevance_table$correction_method[!is.na(result$relevance_table$correction_method)]), "bonferroni")

  cat("\nBonferroni: Selected", result$n_features_selected, "out of", result$n_features_original, "\n")
})

test_that("FWER: Holm correction works (less conservative)", {
  data <- generate_test_data()

  # Holm is uniformly more powerful than Bonferroni
  result <- suppressWarnings({
    select_features(
      data$X, data$y,
      fdr_level = 0.05,
      correction_method = "holm"
    )
  })

  expect_s3_class(result, "tsdiscov_feature_selection")
  expect_gte(result$n_features_selected, 1)
  expect_equal(unique(result$relevance_table$correction_method[!is.na(result$relevance_table$correction_method)]), "holm")

  cat("Holm: Selected", result$n_features_selected, "out of", result$n_features_original, "\n")
})

test_that("FWER: Hochberg correction works", {
  data <- generate_test_data()

  result <- suppressWarnings({
    select_features(
      data$X, data$y,
      fdr_level = 0.05,
      correction_method = "hochberg"
    )
  })

  expect_s3_class(result, "tsdiscov_feature_selection")
  expect_gte(result$n_features_selected, 1)
  expect_equal(unique(result$relevance_table$correction_method[!is.na(result$relevance_table$correction_method)]), "hochberg")

  cat("Hochberg: Selected", result$n_features_selected, "out of", result$n_features_original, "\n")
})

test_that("FWER: Hommel correction works", {
  data <- generate_test_data()

  result <- suppressWarnings({
    select_features(
      data$X, data$y,
      fdr_level = 0.05,
      correction_method = "hommel"
    )
  })

  expect_s3_class(result, "tsdiscov_feature_selection")
  expect_gte(result$n_features_selected, 1)
  expect_equal(unique(result$relevance_table$correction_method[!is.na(result$relevance_table$correction_method)]), "hommel")

  cat("Hommel: Selected", result$n_features_selected, "out of", result$n_features_original, "\n")
})

test_that("FWER methods are more conservative than FDR", {
  data <- generate_test_data()

  # FDR (BY - conservative)
  result_fdr <- suppressWarnings({
    select_features(
      data$X, data$y,
      fdr_level = 0.05,
      correction_method = "fdr",
      hypotheses_independent = FALSE
    )
  })

  # Bonferroni (most conservative FWER)
  result_bonf <- suppressWarnings({
    select_features(
      data$X, data$y,
      fdr_level = 0.05,
      correction_method = "bonferroni"
    )
  })

  # Holm (less conservative FWER)
  result_holm <- suppressWarnings({
    select_features(
      data$X, data$y,
      fdr_level = 0.05,
      correction_method = "holm"
    )
  })

  cat("\nComparison of correction methods:\n")
  cat("  FDR (BY):", result_fdr$n_features_selected, "features\n")
  cat("  Holm:", result_holm$n_features_selected, "features\n")
  cat("  Bonferroni:", result_bonf$n_features_selected, "features\n")

  # FDR should select more features than FWER methods
  expect_gte(result_fdr$n_features_selected, result_bonf$n_features_selected)

  # Holm should select at least as many as Bonferroni
  expect_gte(result_holm$n_features_selected, result_bonf$n_features_selected)
})

# ==============================================================================
# Test Multiclass Classification
# ==============================================================================

test_that("Feature selection works for multiclass classification", {
  set.seed(456)

  # Generate 3 classes
  n_per_class <- 10
  class0 <- lapply(1:n_per_class, function(i) rnorm(100))  # White noise
  class1 <- lapply(1:n_per_class, function(i) arima.sim(list(ar = 0.7), 100))  # AR
  class2 <- lapply(1:n_per_class, function(i) arima.sim(list(ma = 0.7), 100))  # MA

  X <- rbind(
    do.call(rbind, lapply(class0, ts_features)),
    do.call(rbind, lapply(class1, ts_features)),
    do.call(rbind, lapply(class2, ts_features))
  )
  y <- c(rep(0, n_per_class), rep(1, n_per_class), rep(2, n_per_class))

  # Feature selection with multiclass
  result <- suppressWarnings({
    select_features(
      X, y,
      fdr_level = 0.1,
      multiclass = TRUE,
      n_significant = 1
    )
  })

  expect_s3_class(result, "tsdiscov_feature_selection")
  expect_equal(result$ml_task, "classification")
  expect_gte(result$n_features_selected, 1)

  cat("\nMulticlass (3 classes): Selected", result$n_features_selected, "out of", result$n_features_original, "\n")

  # Check for one-vs-all columns
  rel_table <- result$relevance_table
  expect_true(any(grepl("^p_value_", names(rel_table))))
  expect_true(any(grepl("^p_adjusted_", names(rel_table))))
  expect_true("n_significant" %in% names(rel_table))
})

test_that("Multiclass: n_significant parameter works", {
  set.seed(789)

  # Generate 3 classes
  n_per_class <- 10
  class0 <- lapply(1:n_per_class, function(i) rnorm(100))
  class1 <- lapply(1:n_per_class, function(i) cumsum(rnorm(100)))
  class2 <- lapply(1:n_per_class, function(i) 1:100 + rnorm(100, sd = 5))

  X <- rbind(
    do.call(rbind, lapply(class0, ts_features)),
    do.call(rbind, lapply(class1, ts_features)),
    do.call(rbind, lapply(class2, ts_features))
  )
  y <- c(rep(0, n_per_class), rep(1, n_per_class), rep(2, n_per_class))

  # Require significance in at least 2 classes
  result <- suppressWarnings({
    select_features(
      X, y,
      fdr_level = 0.1,
      multiclass = TRUE,
      n_significant = 2
    )
  })

  expect_gte(result$n_features_selected, 1)

  # Check that selected features are significant in >= 2 classes
  selected <- result$relevance_table[result$relevance_table$relevant, ]
  if (nrow(selected) > 0) {
    expect_true(all(selected$n_significant >= 2))
  }

  cat("Multiclass (n_significant=2): Selected", result$n_features_selected, "features\n")
})

test_that("Multiclass with FWER control", {
  set.seed(111)

  n_per_class <- 10
  class0 <- lapply(1:n_per_class, function(i) rnorm(100))
  class1 <- lapply(1:n_per_class, function(i) cumsum(rnorm(100)))
  class2 <- lapply(1:n_per_class, function(i) arima.sim(list(ar = 0.8), 100))

  X <- rbind(
    do.call(rbind, lapply(class0, ts_features)),
    do.call(rbind, lapply(class1, ts_features)),
    do.call(rbind, lapply(class2, ts_features))
  )
  y <- factor(c(rep("WN", n_per_class), rep("RW", n_per_class), rep("AR", n_per_class)))

  # Use Bonferroni correction for multiclass
  result <- suppressWarnings({
    select_features(
      X, y,
      fdr_level = 0.05,
      correction_method = "bonferroni",
      multiclass = TRUE
    )
  })

  expect_s3_class(result, "tsdiscov_feature_selection")
  expect_gte(result$n_features_selected, 1)

  cat("Multiclass + Bonferroni: Selected", result$n_features_selected, "features\n")
})

# ==============================================================================
# Test Regression
# ==============================================================================

test_that("Feature selection works for regression", {
  set.seed(222)

  # Generate time series of different lengths
  # Target: series length (continuous)
  lengths <- c(50, 100, 150, 200, 250)
  ts_list <- lapply(lengths, function(n) {
    lapply(1:5, function(i) rnorm(n))
  })
  ts_all <- unlist(ts_list, recursive = FALSE)

  X <- do.call(rbind, lapply(ts_all, ts_features))
  y <- rep(lengths, each = 5)  # Continuous target

  result <- suppressWarnings({
    select_features(
      X, y,
      fdr_level = 0.1,
      ml_task = "regression"
    )
  })

  expect_s3_class(result, "tsdiscov_feature_selection")
  expect_equal(result$ml_task, "regression")
  expect_gte(result$n_features_selected, 1)

  cat("\nRegression: Selected", result$n_features_selected, "out of", result$n_features_original, "\n")

  # Check that appropriate test was used (Kendall's tau for regression)
  # This would show in how p-values are computed
})

test_that("Regression with FWER control", {
  set.seed(333)

  # Generate time series with different variance levels
  # Target: variance (continuous)
  variances <- c(0.5, 1, 2, 5, 10)
  ts_list <- lapply(variances, function(sd) {
    lapply(1:5, function(i) rnorm(100, sd = sd))
  })
  ts_all <- unlist(ts_list, recursive = FALSE)

  X <- do.call(rbind, lapply(ts_all, ts_features))
  y <- rep(variances, each = 5)

  # Use Holm correction
  result <- suppressWarnings({
    select_features(
      X, y,
      fdr_level = 0.05,
      correction_method = "holm",
      ml_task = "regression"
    )
  })

  expect_s3_class(result, "tsdiscov_feature_selection")
  expect_equal(result$ml_task, "regression")
  expect_gte(result$n_features_selected, 1)

  cat("Regression + Holm: Selected", result$n_features_selected, "features\n")

  # Should select variance-related features
  selected_names <- result$relevance_table$feature[result$relevance_table$relevant]
  variance_features <- grep("std|var|mad|range", selected_names, value = TRUE)
  cat("  Variance-related features:", length(variance_features), "\n")
})

# ==============================================================================
# Error Handling
# ==============================================================================

test_that("Invalid correction method throws error", {
  data <- generate_test_data()

  expect_error(
    select_features(data$X, data$y, correction_method = "invalid"),
    "correction_method must be one of"
  )
})

test_that("All correction methods produce valid output", {
  data <- generate_test_data(n_per_class = 10)

  methods <- c("fdr", "bonferroni", "holm", "hochberg", "hommel")

  for (method in methods) {
    result <- suppressWarnings({
      select_features(
        data$X, data$y,
        fdr_level = 0.05,
        correction_method = method
      )
    })

    expect_s3_class(result, "tsdiscov_feature_selection")
    expect_true(is.numeric(result$n_features_selected))
    expect_true(is.data.frame(result$relevance_table))
    expect_true("relevant" %in% names(result$relevance_table))
  }
})
