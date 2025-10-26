# Test All Feature/Target Type Combinations
# Comprehensive tests for binary/real features × binary/real targets

library(testthat)
library(tsdiscov)

# ==============================================================================
# Test Binary Target + Binary Feature (Fisher's Exact Test)
# ==============================================================================

test_that("Binary target + Binary feature: Fisher's exact test", {
  set.seed(100)

  # Create binary features (0/1)
  # Feature 1: Strongly associated with target
  # Feature 2: Not associated with target
  # Feature 3: Moderately associated

  n <- 50
  y <- c(rep(0, 25), rep(1, 25))

  # Create binary features
  X <- data.frame(
    binary_strong = c(rep(0, 20), rep(1, 5), rep(0, 5), rep(1, 20)),  # Strong
    binary_weak = sample(c(0, 1), n, replace = TRUE),                  # Random
    binary_moderate = c(rep(0, 15), rep(1, 10), rep(0, 10), rep(1, 15)) # Moderate
  )

  result <- suppressWarnings(select_features(X, y, fdr_level = 0.05))

  expect_s3_class(result, "tsdiscov_feature_selection")
  expect_gte(result$n_features_selected, 1)

  # Check that strongly associated feature is selected
  selected <- result$relevance_table$feature[result$relevance_table$relevant]
  expect_true("binary_strong" %in% selected)

  cat("\nBinary target + Binary feature:\n")
  cat("  Selected:", result$n_features_selected, "out of", result$n_features_original, "\n")
  cat("  Selected features:", paste(selected, collapse = ", "), "\n")
})

# ==============================================================================
# Test Binary Target + Real Feature (Mann-Whitney U Test)
# ==============================================================================

test_that("Binary target + Real feature: Mann-Whitney U test", {
  set.seed(101)

  n_per_class <- 25
  y <- c(rep(0, n_per_class), rep(1, n_per_class))

  # Create real-valued features with different discriminative power
  X <- data.frame(
    real_strong = c(rnorm(n_per_class, mean = 0, sd = 1),
                    rnorm(n_per_class, mean = 2, sd = 1)),    # Strong separation
    real_weak = rnorm(2 * n_per_class, mean = 0, sd = 1),    # No separation
    real_moderate = c(rnorm(n_per_class, mean = 0, sd = 1.5),
                      rnorm(n_per_class, mean = 1, sd = 1.5))  # Moderate
  )

  result <- suppressWarnings(select_features(X, y, fdr_level = 0.05))

  expect_s3_class(result, "tsdiscov_feature_selection")
  expect_gte(result$n_features_selected, 1)

  # Check that strongly discriminative feature is selected
  selected <- result$relevance_table$feature[result$relevance_table$relevant]
  expect_true("real_strong" %in% selected)

  cat("\nBinary target + Real feature:\n")
  cat("  Selected:", result$n_features_selected, "out of", result$n_features_original, "\n")
  cat("  P-values:\n")
  print(result$relevance_table[, c("feature", "type", "p_value")])
})

test_that("Binary target + Real feature: KS test option", {
  set.seed(102)

  n_per_class <- 25
  y <- c(rep(0, n_per_class), rep(1, n_per_class))

  X <- data.frame(
    feature1 = c(rnorm(n_per_class, 0, 1), rnorm(n_per_class, 1.5, 1))
  )

  # Use KS test instead of Mann-Whitney
  result <- suppressWarnings(
    select_features(X, y, fdr_level = 0.05, test_binary_real = "ks")
  )

  expect_s3_class(result, "tsdiscov_feature_selection")
  expect_gte(result$n_features_selected, 0)

  cat("\nBinary target + Real feature (KS test):\n")
  cat("  Selected:", result$n_features_selected, "\n")
})

# ==============================================================================
# Test Real Target + Binary Feature (KS Test)
# ==============================================================================

test_that("Real target + Binary feature: KS test", {
  set.seed(103)

  n <- 50
  y <- rnorm(n, mean = 10, sd = 3)  # Real-valued target

  # Create binary features
  # Feature 1: Values differ based on target magnitude
  X <- data.frame(
    binary_associated = ifelse(y > median(y), 1, 0),  # Associated
    binary_random = sample(c(0, 1), n, replace = TRUE)  # Random
  )

  result <- suppressWarnings(
    select_features(X, y, fdr_level = 0.05, ml_task = "regression")
  )

  expect_s3_class(result, "tsdiscov_feature_selection")
  expect_equal(result$ml_task, "regression")

  cat("\nReal target + Binary feature:\n")
  cat("  Selected:", result$n_features_selected, "out of", result$n_features_original, "\n")
  if (result$n_features_selected > 0) {
    selected <- result$relevance_table$feature[result$relevance_table$relevant]
    cat("  Selected features:", paste(selected, collapse = ", "), "\n")
  }
})

# ==============================================================================
# Test Real Target + Real Feature (Kendall's Tau)
# ==============================================================================

test_that("Real target + Real feature: Kendall's tau", {
  set.seed(104)

  n <- 50
  y <- rnorm(n, mean = 100, sd = 20)  # Real-valued target

  # Create real features with different correlations
  X <- data.frame(
    strong_positive = y + rnorm(n, 0, 5),      # Strong positive correlation
    weak = rnorm(n, 50, 10),                    # No correlation
    moderate_negative = -0.5 * y + rnorm(n, 0, 15)  # Moderate negative
  )

  result <- suppressWarnings(
    select_features(X, y, fdr_level = 0.05, ml_task = "regression")
  )

  expect_s3_class(result, "tsdiscov_feature_selection")
  expect_equal(result$ml_task, "regression")
  expect_gte(result$n_features_selected, 1)

  # Strongly correlated features should be selected
  selected <- result$relevance_table$feature[result$relevance_table$relevant]
  expect_true("strong_positive" %in% selected)

  cat("\nReal target + Real feature:\n")
  cat("  Selected:", result$n_features_selected, "out of", result$n_features_original, "\n")
  cat("  P-values:\n")
  print(result$relevance_table[, c("feature", "type", "p_value")])
})

# ==============================================================================
# Test Mixed Feature Types with Binary Target
# ==============================================================================

test_that("Mixed feature types with binary target", {
  set.seed(105)

  n_per_class <- 30
  y <- c(rep(0, n_per_class), rep(1, n_per_class))

  # Mix of binary and real features
  X <- data.frame(
    binary1 = c(rep(0, 25), rep(1, 5), rep(0, 5), rep(1, 25)),
    real1 = c(rnorm(n_per_class, 0, 1), rnorm(n_per_class, 2, 1)),
    binary2 = sample(c(0, 1), 2 * n_per_class, replace = TRUE),
    real2 = rnorm(2 * n_per_class, 10, 3),
    real3 = c(rnorm(n_per_class, 5, 2), rnorm(n_per_class, 8, 2))
  )

  result <- suppressWarnings(select_features(X, y, fdr_level = 0.05))

  expect_s3_class(result, "tsdiscov_feature_selection")
  expect_gte(result$n_features_selected, 2)

  # Check feature types in relevance table
  expect_true("type" %in% names(result$relevance_table))
  types <- unique(result$relevance_table$type)
  expect_true(any(types %in% c("binary", "real")))

  cat("\nMixed feature types (binary target):\n")
  cat("  Total features:", result$n_features_original, "\n")
  cat("  Selected:", result$n_features_selected, "\n")
  cat("  Feature types:\n")
  print(table(result$relevance_table$type))
  cat("  Selected by type:\n")
  selected_table <- result$relevance_table[result$relevance_table$relevant, ]
  if (nrow(selected_table) > 0) {
    print(table(selected_table$type))
  }
})

# ==============================================================================
# Test Mixed Feature Types with Real Target
# ==============================================================================

test_that("Mixed feature types with real target", {
  set.seed(106)

  n <- 60
  y <- seq(50, 150, length.out = n) + rnorm(n, 0, 10)  # Continuous target with trend

  # Mix of binary and real features
  X <- data.frame(
    binary1 = ifelse(y > median(y), 1, 0),
    real1 = 0.8 * y + rnorm(n, 0, 10),
    binary2 = sample(c(0, 1), n, replace = TRUE),
    real2 = rnorm(n, 100, 20),
    real3 = -0.6 * y + rnorm(n, 150, 15)
  )

  result <- suppressWarnings(
    select_features(X, y, fdr_level = 0.05, ml_task = "regression")
  )

  expect_s3_class(result, "tsdiscov_feature_selection")
  expect_equal(result$ml_task, "regression")
  expect_gte(result$n_features_selected, 1)

  cat("\nMixed feature types (real target):\n")
  cat("  Selected:", result$n_features_selected, "out of", result$n_features_original, "\n")
  cat("  P-values:\n")
  print(result$relevance_table[, c("feature", "type", "p_value")])
})

# ==============================================================================
# Test Constant Features (Edge Case)
# ==============================================================================

test_that("Constant features are never selected", {
  set.seed(107)

  n <- 50
  y <- c(rep(0, 25), rep(1, 25))

  X <- data.frame(
    constant = rep(5, n),                              # Constant
    good_feature = c(rnorm(25, 0, 1), rnorm(25, 2, 1)), # Discriminative
    nearly_constant = c(rep(1, 49), 2)                  # Nearly constant
  )

  result <- suppressWarnings(select_features(X, y, fdr_level = 0.05))

  # Constant feature should never be selected
  expect_false("constant" %in% result$relevance_table$feature[result$relevance_table$relevant])

  # Good feature should be selected
  expect_true("good_feature" %in% result$relevance_table$feature[result$relevance_table$relevant])

  cat("\nConstant features:\n")
  cat("  Feature types:\n")
  print(result$relevance_table[, c("feature", "type", "relevant")])
})

# ==============================================================================
# Test All Combinations with FWER Control
# ==============================================================================

test_that("All type combinations work with Bonferroni correction", {
  set.seed(108)

  test_cases <- list(
    list(
      name = "Binary target + Binary feature",
      X = data.frame(f = c(rep(0, 25), rep(1, 25))),
      y = c(rep(0, 25), rep(1, 25)),
      ml_task = "auto"
    ),
    list(
      name = "Binary target + Real feature",
      X = data.frame(f = c(rnorm(25, 0), rnorm(25, 2))),
      y = c(rep(0, 25), rep(1, 25)),
      ml_task = "auto"
    ),
    list(
      name = "Real target + Binary feature",
      X = data.frame(f = rep(c(0, 1), each = 25)),
      y = c(rnorm(25, 10), rnorm(25, 20)),
      ml_task = "regression"
    ),
    list(
      name = "Real target + Real feature",
      X = data.frame(f = seq(1, 50)),
      y = seq(10, 100, length.out = 50) + rnorm(50, 0, 5),
      ml_task = "regression"
    )
  )

  for (tc in test_cases) {
    result <- suppressWarnings(
      select_features(
        tc$X, tc$y,
        fdr_level = 0.05,
        correction_method = "bonferroni",
        ml_task = tc$ml_task
      )
    )

    expect_s3_class(result, "tsdiscov_feature_selection")
    cat(sprintf("\n%s: ", tc$name))
    cat(result$n_features_selected, "selected,",
        "p-value =", format(result$relevance_table$p_value[1], digits = 3), "\n")
  }
})

# ==============================================================================
# Summary Test: Type Detection
# ==============================================================================

test_that("Feature type classification works correctly", {
  set.seed(109)

  n <- 50
  X <- data.frame(
    constant_feature = rep(1, n),
    binary_feature = sample(c(0, 1), n, replace = TRUE),
    three_values = sample(c(1, 2, 3), n, replace = TRUE),
    real_feature = rnorm(n)
  )
  y <- c(rep(0, 25), rep(1, 25))

  result <- suppressWarnings(select_features(X, y, fdr_level = 0.1))

  # Check type classification
  rel_table <- result$relevance_table
  expect_equal(rel_table$type[rel_table$feature == "constant_feature"], "constant")
  expect_equal(rel_table$type[rel_table$feature == "binary_feature"], "binary")
  expect_equal(rel_table$type[rel_table$feature == "three_values"], "real")
  expect_equal(rel_table$type[rel_table$feature == "real_feature"], "real")

  cat("\nFeature type detection:\n")
  print(rel_table[, c("feature", "type")])
})
