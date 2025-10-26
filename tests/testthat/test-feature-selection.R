# Tests for Feature Selection
# Implements FRESH algorithm (FeatuRE Selection with Hypothesis testing)

library(testthat)
library(tsdiscov)

# =============================================================================
# Helper Functions
# =============================================================================

generate_test_data <- function(n = 100, n_relevant = 3, n_noise = 5) {
  set.seed(123)

  # Binary classification
  y <- c(rep(0, n/2), rep(1, n/2))

  # Relevant features
  X <- data.frame(
    relevant1 = c(rnorm(n/2, 0, 1), rnorm(n/2, 2, 1)),  # Shifted mean
    relevant2 = c(rnorm(n/2, -1, 1), rnorm(n/2, 1, 1)), # Shifted mean
    relevant3 = c(rnorm(n/2, 0, 0.5), rnorm(n/2, 0, 2)) # Different variance
  )

  # Noise features
  for (i in seq_len(n_noise)) {
    X[[paste0("noise", i)]] <- rnorm(n)
  }

  list(X = X, y = y)
}

# =============================================================================
# Tests for Feature Classification
# =============================================================================

test_that("classify_feature identifies constant features", {
  expect_equal(classify_feature(rep(5, 10)), "constant")
  expect_equal(classify_feature(rep(0, 100)), "constant")
})

test_that("classify_feature identifies binary features", {
  expect_equal(classify_feature(c(0, 1, 0, 1, 0)), "binary")
  expect_equal(classify_feature(c(TRUE, FALSE, TRUE)), "binary")
  expect_equal(classify_feature(c(10, 20, 10, 20)), "binary")
})

test_that("classify_feature identifies real features", {
  expect_equal(classify_feature(rnorm(100)), "real")
  expect_equal(classify_feature(1:10), "real")
  expect_equal(classify_feature(c(1, 2, 3, 4, 5)), "real")
})

test_that("classify_feature handles NA values", {
  expect_equal(classify_feature(c(1, 2, NA, NA)), "binary")
  expect_equal(classify_feature(c(1, NA, NA)), "constant")
  expect_equal(classify_feature(c(1, 2, 3, NA)), "real")
})

# =============================================================================
# Tests for ML Task Inference
# =============================================================================

test_that("infer_ml_task identifies classification", {
  expect_equal(infer_ml_task(factor(c("A", "B", "A"))), "classification")
  expect_equal(infer_ml_task(c(TRUE, FALSE, TRUE)), "classification")
  expect_equal(infer_ml_task(c(0L, 1L, 0L, 1L)), "classification")
  expect_equal(infer_ml_task(c("cat", "dog", "cat")), "classification")
})

test_that("infer_ml_task identifies regression", {
  expect_equal(infer_ml_task(rnorm(100)), "regression")
  expect_equal(infer_ml_task(seq(0, 10, 0.1)), "regression")
  expect_equal(infer_ml_task(runif(50)), "regression")
})

# =============================================================================
# Tests for Statistical Tests
# =============================================================================

test_that("test_binary_binary works correctly", {
  # Perfect separation
  feature <- c(rep(0, 50), rep(1, 50))
  target <- c(rep(0, 50), rep(1, 50))
  p <- test_binary_binary(feature, target)
  expect_true(is.numeric(p))
  expect_true(p < 0.01)  # Should be highly significant

  # No association
  set.seed(123)
  feature <- sample(c(0, 1), 100, replace = TRUE)
  target <- sample(c(0, 1), 100, replace = TRUE)
  p <- test_binary_binary(feature, target)
  expect_true(is.numeric(p))
  expect_true(p > 0.05)  # Should not be significant
})

test_that("test_binary_real works with Mann-Whitney U", {
  # Clear difference
  feature <- c(rnorm(50, 0), rnorm(50, 2))
  target <- c(rep(0, 50), rep(1, 50))
  p <- test_binary_real(feature, target, test = "mann")
  expect_true(is.numeric(p))
  expect_true(p < 0.01)

  # No difference
  set.seed(123)
  feature <- rnorm(100)
  target <- sample(c(0, 1), 100, replace = TRUE)
  p <- test_binary_real(feature, target, test = "mann")
  expect_true(is.numeric(p))
  expect_true(p > 0.05)
})

test_that("test_binary_real works with Kolmogorov-Smirnov", {
  feature <- c(rnorm(50, 0), rnorm(50, 2))
  target <- c(rep(0, 50), rep(1, 50))
  p <- test_binary_real(feature, target, test = "ks")
  expect_true(is.numeric(p))
  expect_true(p < 0.01)
})

test_that("test_real_binary works correctly", {
  feature <- c(rep(0, 50), rep(1, 50))
  target <- c(rnorm(50, 0), rnorm(50, 2))
  p <- test_real_binary(feature, target)
  expect_true(is.numeric(p))
  expect_true(p < 0.01)
})

test_that("test_real_real works correctly", {
  # Strong correlation
  set.seed(123)
  feature <- 1:100
  target <- 2 * feature + rnorm(100, sd = 5)
  p <- test_real_real(feature, target)
  expect_true(is.numeric(p))
  expect_true(p < 0.01)

  # No correlation
  set.seed(456)
  feature <- rnorm(100)
  target <- rnorm(100)
  p <- test_real_real(feature, target)
  expect_true(is.numeric(p))
  expect_true(p > 0.05)
})

# =============================================================================
# Tests for calculate_relevance_table
# =============================================================================

test_that("calculate_relevance_table works for binary classification", {
  data <- generate_test_data()

  relevance <- calculate_relevance_table(
    data$X, data$y,
    fdr_level = 0.05,
    ml_task = "classification"
  )

  expect_true(is.data.frame(relevance))
  expect_equal(nrow(relevance), ncol(data$X))
  expect_true(all(c("feature", "type", "p_value", "p_adjusted", "relevant") %in% names(relevance)))

  # Relevant features should have low p-values
  relevant_features <- relevance$feature[relevance$relevant]
  expect_true("relevant1" %in% relevant_features)
  expect_true("relevant2" %in% relevant_features)
})

test_that("calculate_relevance_table works for regression", {
  set.seed(123)
  X <- data.frame(
    relevant = seq(0, 10, length.out = 100),
    noise = rnorm(100)
  )
  y <- 2 * X$relevant + rnorm(100, sd = 1)

  relevance <- calculate_relevance_table(
    X, y,
    fdr_level = 0.05,
    ml_task = "regression"
  )

  expect_true("relevant" %in% relevance$feature[relevance$relevant])
  expect_false("noise" %in% relevance$feature[relevance$relevant])
})

test_that("calculate_relevance_table handles constant features", {
  X <- data.frame(
    constant = rep(5, 100),
    variable = rnorm(100)
  )
  y <- sample(c(0, 1), 100, replace = TRUE)

  relevance <- calculate_relevance_table(X, y, ml_task = "classification")

  const_row <- relevance[relevance$feature == "constant", ]
  expect_equal(const_row$type, "constant")
  expect_false(const_row$relevant)
})

test_that("calculate_relevance_table uses FDR correction", {
  set.seed(123)
  # Many noise features
  X <- data.frame(replicate(50, rnorm(100)))
  y <- sample(c(0, 1), 100, replace = TRUE)

  relevance <- calculate_relevance_table(X, y, fdr_level = 0.05)

  # With FDR, should select few features (expected ~2.5)
  expect_true(sum(relevance$relevant) < 10)

  # p_adjusted should be >= p_value (correction makes them more conservative)
  valid <- !is.na(relevance$p_value)
  expect_true(all(relevance$p_adjusted[valid] >= relevance$p_value[valid]))
})

# =============================================================================
# Tests for select_features
# =============================================================================

test_that("select_features basic functionality", {
  data <- generate_test_data()

  result <- select_features(data$X, data$y, fdr_level = 0.05)

  expect_s3_class(result, "tsdiscov_feature_selection")
  expect_true(is.data.frame(result$X_selected))
  expect_true(is.data.frame(result$relevance_table))
  expect_true(is.numeric(result$n_features_original))
  expect_true(is.numeric(result$n_features_selected))
})

test_that("select_features filters correctly", {
  data <- generate_test_data(n = 100, n_relevant = 3, n_noise = 10)

  result <- select_features(data$X, data$y, fdr_level = 0.05)

  # Should select relevant features
  expect_true("relevant1" %in% colnames(result$X_selected))
  expect_true("relevant2" %in% colnames(result$X_selected))

  # Should have fewer features than original
  expect_true(result$n_features_selected < result$n_features_original)
  expect_equal(ncol(result$X_selected), result$n_features_selected)
})

test_that("select_features works with matrix input", {
  data <- generate_test_data()
  X_matrix <- as.matrix(data$X)

  result <- select_features(X_matrix, data$y)

  expect_s3_class(result, "tsdiscov_feature_selection")
  expect_true(is.data.frame(result$X_selected))
})

test_that("select_features works with factor target", {
  data <- generate_test_data()
  y_factor <- factor(data$y, labels = c("class0", "class1"))

  result <- select_features(data$X, y_factor)

  expect_s3_class(result, "tsdiscov_feature_selection")
  expect_equal(result$ml_task, "classification")
})

test_that("select_features auto-detects classification", {
  data <- generate_test_data()

  result <- select_features(data$X, data$y, ml_task = "auto")

  expect_equal(result$ml_task, "classification")
})

test_that("select_features auto-detects regression", {
  set.seed(123)
  X <- data.frame(feature = rnorm(100))
  y <- rnorm(100)

  result <- select_features(X, y, ml_task = "auto")

  expect_equal(result$ml_task, "regression")
})

# =============================================================================
# Tests for Input Validation
# =============================================================================

test_that("select_features validates inputs", {
  data <- generate_test_data()

  # Not a data.frame/matrix
  expect_error(select_features(list(a = 1, b = 2), data$y))

  # Mismatched dimensions
  expect_error(select_features(data$X, data$y[1:50]))

  # Too few observations
  expect_error(select_features(data$X[1, , drop = FALSE], data$y[1]))

  # Only one class
  expect_error(select_features(data$X, rep(0, 100)))

  # Invalid fdr_level
  expect_error(select_features(data$X, data$y, fdr_level = -0.1))
  expect_error(select_features(data$X, data$y, fdr_level = 1.5))
})

# =============================================================================
# Tests for Multiclass Classification
# =============================================================================

test_that("multiclass classification works", {
  set.seed(123)
  N <- 30

  # 3-class problem
  X <- data.frame(
    f1 = c(rnorm(N, 0), rnorm(N, 2), rnorm(N, 0)),   # Discriminates class 1
    f2 = c(rnorm(N, 0), rnorm(N, 0), rnorm(N, 2)),   # Discriminates class 2
    f3 = c(rnorm(N, 1), rnorm(N, 1), rnorm(N, 1))    # Noise
  )
  y <- rep(0:2, each = N)

  result <- select_features(
    X, y,
    multiclass = TRUE,
    n_significant = 1,
    fdr_level = 0.05
  )

  expect_s3_class(result, "tsdiscov_feature_selection")

  # Should select f1 and f2 (each discriminates at least 1 class)
  expect_true("f1" %in% colnames(result$X_selected))
  expect_true("f2" %in% colnames(result$X_selected))

  # Relevance table should have per-class columns
  expect_true(any(grepl("^p_value_", names(result$relevance_table))))
  expect_true(any(grepl("^relevant_", names(result$relevance_table))))
  expect_true("n_significant" %in% names(result$relevance_table))
})

test_that("multiclass n_significant parameter works", {
  # Test multiclass feature selection with controlled feature relevance
  # Key insight: Start with zeros, then replace specific class values for perfect separation
  set.seed(456)
  N <- 100  # Samples per class

  # Create 3-class target
  y <- rep(0:2, each = N)

  # Create features with controlled n_significant values
  # Start with mostly zeros
  base <- rep(0, N * 3)

  # f_none: All zeros (no separation) -> n_significant = 0
  f_none <- base

  # f_one: Only class 0 is different -> n_significant should vary
  # (might be 1 or 3 depending on one-vs-all dynamics)
  f_one <- base
  f_one[y == 0] <- runif(N, 2, 3)  # Class 0: [2,3], others: 0

  # f_two: Classes 0 and 1 are different (in opposite directions) -> n_significant = 2
  # Class 2 distribution (all zeros) is contained in the range of "all others"
  # so class 2 vs rest will NOT be significant
  f_two <- base
  f_two[y == 0] <- runif(N, 2, 3)   # Class 0: [2,3]
  f_two[y == 1] <- runif(N, -2, -1) # Class 1: [-2,-1]
  # Class 2: 0 (contained in the range of classes 0 and 1)

  X <- data.frame(f_none = f_none, f_one = f_one, f_two = f_two)

  # Test with n_significant = 1: should select f_one and f_two
  result1 <- select_features(X, y, multiclass = TRUE, n_significant = 1, fdr_level = 0.05)
  expect_true("f_one" %in% colnames(result1$X_selected))
  expect_true("f_two" %in% colnames(result1$X_selected))
  expect_false("f_none" %in% colnames(result1$X_selected))

  # Test with n_significant = 2: should select only f_two
  result2 <- select_features(X, y, multiclass = TRUE, n_significant = 2, fdr_level = 0.05)
  expect_true("f_two" %in% colnames(result2$X_selected),
              info = "f_two should be selected (n_significant=2)")

  # Verify n_significant counts in relevance table
  rel_table <- result2$relevance_table
  expect_equal(rel_table$n_significant[rel_table$feature == "f_two"], 2,
               info = "f_two should be significant for exactly 2 classes")
})

# =============================================================================
# Tests for Print/Summary Methods
# =============================================================================

test_that("print method works", {
  data <- generate_test_data()
  result <- select_features(data$X, data$y)

  expect_output(print(result), "Feature Selection Results")
  expect_output(print(result), "Original features:")
  expect_output(print(result), "Selected features:")
})

test_that("summary method works", {
  data <- generate_test_data()
  result <- select_features(data$X, data$y)

  expect_output(summary(result), "Feature Selection Summary")
  expect_output(summary(result), "Features by type:")
  expect_output(summary(result), "P-value distribution:")
})

# =============================================================================
# Integration Tests
# =============================================================================

test_that("Integration with ts_features", {
  set.seed(123)

  # Generate different types of time series
  ts_normal <- lapply(1:20, function(i) rnorm(100))
  ts_ar <- lapply(1:20, function(i) arima.sim(list(ar = 0.8), 100))

  # Extract features (using a subset for speed)
  features_normal <- do.call(rbind, lapply(ts_normal, function(x) {
    ts_features(x, features = c("stats", "acf", "trend"))
  }))

  features_ar <- do.call(rbind, lapply(ts_ar, function(x) {
    ts_features(x, features = c("stats", "acf", "trend"))
  }))

  X <- rbind(features_normal, features_ar)
  y <- c(rep(0, 20), rep(1, 20))

  # Select features
  result <- select_features(X, y, fdr_level = 0.10)

  expect_s3_class(result, "tsdiscov_feature_selection")
  expect_true(result$n_features_selected > 0)
  expect_true(result$n_features_selected < result$n_features_original)

  # ACF features should be selected (AR process has different autocorrelation)
  relevant_features <- colnames(result$X_selected)
  acf_selected <- any(grepl("acf", relevant_features))
  expect_true(acf_selected)
})

test_that("Feature selection with all 210 features", {
  skip_if_not(interactive(), "Skipping slow test in non-interactive mode")

  set.seed(456)

  # Generate time series
  ts1 <- lapply(1:10, function(i) rnorm(100))
  ts2 <- lapply(1:10, function(i) arima.sim(list(ar = 0.7), 100))

  # Extract all features
  features1 <- do.call(rbind, lapply(ts1, ts_features))
  features2 <- do.call(rbind, lapply(ts2, ts_features))

  X <- rbind(features1, features2)
  y <- c(rep(0, 10), rep(1, 10))

  # Select features
  result <- select_features(X, y, fdr_level = 0.10)

  expect_true(result$n_features_selected < 213)
  expect_true(result$n_features_selected > 0)
})

# =============================================================================
# Edge Cases
# =============================================================================

test_that("Handles features with NA values", {
  X <- data.frame(
    good = rnorm(100),
    has_na = c(rnorm(50), rep(NA, 50))
  )
  y <- sample(c(0, 1), 100, replace = TRUE)

  expect_warning(
    result <- select_features(X, y),
    "contains NA"
  )

  expect_s3_class(result, "tsdiscov_feature_selection")
})

test_that("Handles all constant features", {
  X <- data.frame(
    const1 = rep(1, 100),
    const2 = rep(5, 100)
  )
  y <- sample(c(0, 1), 100, replace = TRUE)

  result <- select_features(X, y)

  expect_equal(result$n_features_selected, 0)
  expect_equal(nrow(result$X_selected), 100)
  expect_equal(ncol(result$X_selected), 0)
})

test_that("Handles mixed feature types", {
  set.seed(123)
  X <- data.frame(
    binary1 = sample(c(0, 1), 100, replace = TRUE),
    binary2 = sample(c(0, 1), 100, replace = TRUE),
    real1 = rnorm(100),
    real2 = rnorm(100),
    constant = rep(5, 100)
  )
  y <- sample(c(0, 1), 100, replace = TRUE)

  result <- select_features(X, y)

  expect_s3_class(result, "tsdiscov_feature_selection")

  # Check feature types in relevance table
  types <- unique(result$relevance_table$type)
  expect_true("binary" %in% types)
  expect_true("real" %in% types)
  expect_true("constant" %in% types)
})

test_that("Different FDR levels affect selection", {
  data <- generate_test_data(n = 100, n_relevant = 2, n_noise = 20)

  result_strict <- select_features(data$X, data$y, fdr_level = 0.01)
  result_lenient <- select_features(data$X, data$y, fdr_level = 0.10)

  # More lenient FDR should select more (or equal) features
  expect_true(result_lenient$n_features_selected >= result_strict$n_features_selected)
})

test_that("hypotheses_independent parameter works", {
  data <- generate_test_data()

  result_dep <- select_features(data$X, data$y, hypotheses_independent = FALSE)
  result_indep <- select_features(data$X, data$y, hypotheses_independent = TRUE)

  # Independent assumption should be less conservative (select more features)
  expect_true(result_indep$n_features_selected >= result_dep$n_features_selected)
})
