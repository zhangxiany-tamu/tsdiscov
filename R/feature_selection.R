#' Feature Selection Using Hypothesis Testing
#'
#' Implements the FRESH (FeatuRE Selection with Hypothesis testing) algorithm
#' to automatically select relevant features based on statistical significance
#' tests with False Discovery Rate (FDR) correction.
#'
#' @param X Feature matrix (data.frame or matrix with named columns)
#' @param y Target vector (numeric, factor, or logical)
#' @param fdr_level False discovery rate threshold (default 0.05). Also used as
#'   FWER threshold when correction_method is "bonferroni", "holm", etc.
#' @param hypotheses_independent Assume independent features? If FALSE (default),
#'   uses more conservative Benjamini-Yekutieli correction (only for FDR methods)
#' @param correction_method Multiple testing correction method. Options:
#'   \itemize{
#'     \item "fdr": False Discovery Rate using BH or BY (default)
#'     \item "bonferroni": FWER control, very conservative, robust to any dependence
#'     \item "holm": FWER control, less conservative than Bonferroni
#'     \item "hochberg": FWER control, assumes non-negative dependence
#'     \item "hommel": FWER control, more powerful than Hochberg
#'   }
#' @param ml_task "auto" (default), "classification", or "regression"
#' @param multiclass Is this a multiclass classification problem? (default FALSE)
#' @param n_significant Minimum number of classes a feature must be significant
#'   for in multiclass problems (default 1)
#' @param test_binary_real Test for binary target + real feature: "mann" (default)
#'   for Mann-Whitney U or "ks" for Kolmogorov-Smirnov
#'
#' @return A list of class 'tsdiscov_feature_selection' with components:
#' \itemize{
#'   \item X_selected: Reduced feature matrix with only relevant features
#'   \item relevance_table: data.frame with detailed results for each feature
#'   \item n_features_original: Original number of features
#'   \item n_features_selected: Number of selected features
#'   \item fdr_level: FDR threshold used
#'   \item ml_task: Machine learning task (classification or regression)
#' }
#'
#' @details
#' The FRESH algorithm tests each feature independently for association with
#' the target variable using appropriate statistical tests:
#'
#' \itemize{
#'   \item Binary target + Binary feature: Fisher's exact test
#'   \item Binary target + Real feature: Mann-Whitney U test (default) or
#'     Kolmogorov-Smirnov test
#'   \item Real target + Binary feature: Kolmogorov-Smirnov test
#'   \item Real target + Real feature: Kendall's tau correlation test
#' }
#'
#' Multiple testing correction is applied using the Benjamini-Hochberg (BH)
#' procedure if hypotheses_independent=TRUE, or the more conservative
#' Benjamini-Yekutieli (BY) procedure if hypotheses_independent=FALSE (recommended
#' for correlated features like time series features).
#'
#' @references
#' Christ, M., Kempa-Liehr, A.W. and Feindt, M. (2016).
#' Distributed and parallel time series feature extraction for industrial big data applications.
#' ArXiv e-prints: 1610.07717
#'
#' Benjamini, Y. and Hochberg, Y. (1995).
#' Controlling the false discovery rate: a practical and powerful approach to multiple testing.
#' Journal of the Royal Statistical Society. Series B, 57(1), 289-300.
#'
#' @importFrom stats fisher.test wilcox.test ks.test cor.test p.adjust
#' @export
#' @examples
#' \dontrun{
#' # Extract features from time series
#' library(tsdiscov)
#' set.seed(123)
#'
#' # Generate time series with different characteristics
#' ts_normal <- replicate(50, rnorm(100), simplify = FALSE)
#' ts_ar <- replicate(50, arima.sim(list(ar = 0.8), 100), simplify = FALSE)
#'
#' # Extract all features
#' features_normal <- do.call(rbind, lapply(ts_normal, ts_features))
#' features_ar <- do.call(rbind, lapply(ts_ar, ts_features))
#'
#' X <- rbind(features_normal, features_ar)
#' y <- c(rep(0, 50), rep(1, 50))
#'
#' # Select relevant features
#' result <- select_features(X, y, fdr_level = 0.05)
#'
#' # Examine results
#' print(result)
#' head(result$relevance_table)
#'
#' # Use selected features for modeling
#' X_selected <- result$X_selected
#' }
select_features <- function(X, y,
                           fdr_level = 0.05,
                           hypotheses_independent = FALSE,
                           correction_method = "fdr",
                           ml_task = "auto",
                           multiclass = FALSE,
                           n_significant = 1,
                           test_binary_real = "mann") {

  # Input validation
  if (!is.data.frame(X) && !is.matrix(X)) {
    stop("X must be a data.frame or matrix")
  }

  if (!is.numeric(y) && !is.factor(y) && !is.logical(y)) {
    stop("y must be numeric, factor, or logical")
  }

  if (nrow(X) != length(y)) {
    stop("X and y must have the same number of observations")
  }

  if (nrow(X) < 2) {
    stop("X must have at least 2 observations")
  }

  if (length(unique(y)) < 2) {
    stop("y must have at least 2 unique values")
  }

  if (fdr_level <= 0 || fdr_level >= 1) {
    stop("fdr_level must be between 0 and 1")
  }

  # Validate correction method
  valid_methods <- c("fdr", "bonferroni", "holm", "hochberg", "hommel")
  if (!correction_method %in% valid_methods) {
    stop("correction_method must be one of: ", paste(valid_methods, collapse = ", "))
  }

  # Convert matrix to data.frame if needed
  if (is.matrix(X)) {
    X <- as.data.frame(X)
    if (is.null(colnames(X))) {
      colnames(X) <- paste0("feature_", seq_len(ncol(X)))
    }
  }

  # Ensure X has column names
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("feature_", seq_len(ncol(X)))
  }

  # Calculate relevance table
  if (multiclass) {
    relevance <- calculate_relevance_multiclass(
      X, y,
      fdr_level = fdr_level,
      hypotheses_independent = hypotheses_independent,
      correction_method = correction_method,
      n_significant = n_significant,
      test_binary_real = test_binary_real
    )
  } else {
    relevance <- calculate_relevance_table(
      X, y,
      fdr_level = fdr_level,
      hypotheses_independent = hypotheses_independent,
      correction_method = correction_method,
      ml_task = ml_task,
      test_binary_real = test_binary_real
    )
  }

  # Filter to relevant features
  relevant_features <- relevance$feature[relevance$relevant]
  X_selected <- X[, relevant_features, drop = FALSE]

  # Create result object
  result <- structure(
    list(
      X_selected = X_selected,
      relevance_table = relevance,
      n_features_original = ncol(X),
      n_features_selected = length(relevant_features),
      fdr_level = fdr_level,
      ml_task = if (multiclass) "classification" else
                if (ml_task == "auto") infer_ml_task(y) else ml_task
    ),
    class = "tsdiscov_feature_selection"
  )

  return(result)
}


#' Calculate Relevance Table for Features
#'
#' @param X Feature matrix
#' @param y Target vector
#' @param fdr_level FDR threshold
#' @param hypotheses_independent Independence assumption
#' @param correction_method Multiple testing correction method
#' @param ml_task Machine learning task
#' @param test_binary_real Test for binary target + real feature
#'
#' @return data.frame with relevance information
#' @keywords internal
calculate_relevance_table <- function(X, y,
                                     fdr_level = 0.05,
                                     hypotheses_independent = FALSE,
                                     correction_method = "fdr",
                                     ml_task = "auto",
                                     test_binary_real = "mann") {

  # Infer ML task if needed
  if (ml_task == "auto") {
    ml_task <- infer_ml_task(y)
  }

  # Classify each feature
  feature_types <- vapply(X, classify_feature, character(1))

  # Initialize results
  results <- data.frame(
    feature = colnames(X),
    type = feature_types,
    p_value = NA_real_,
    stringsAsFactors = FALSE
  )

  # Test each feature
  for (i in seq_along(X)) {
    feature <- X[[i]]
    ftype <- feature_types[i]

    # Skip constant features
    if (ftype == "constant") {
      next
    }

    # Check for NAs
    if (any(is.na(feature)) || any(is.na(y))) {
      warning(sprintf("Feature '%s' or target contains NA values, skipping", results$feature[i]))
      next
    }

    # Convert feature to numeric if it's not
    if (!is.numeric(feature)) {
      feature <- as.numeric(feature)
    }

    # Select and run appropriate test
    tryCatch({
      if (ml_task == "classification") {
        # Convert y to binary if needed
        y_levels <- unique(y)
        if (length(y_levels) > 2) {
          stop("Use multiclass=TRUE for multiclass classification")
        }
        y_binary <- as.numeric(factor(y, levels = y_levels)) - 1

        if (ftype == "binary") {
          results$p_value[i] <- test_binary_binary(feature, y_binary)
        } else {
          results$p_value[i] <- test_binary_real(feature, y_binary, test = test_binary_real)
        }
      } else {  # regression
        if (ftype == "binary") {
          results$p_value[i] <- test_real_binary(feature, y)
        } else {
          results$p_value[i] <- test_real_real(feature, y)
        }
      }
    }, error = function(e) {
      warning(sprintf("Error testing feature '%s': %s", results$feature[i], e$message))
    })
  }

  # Apply multiple testing correction
  valid_pvalues <- !is.na(results$p_value)

  if (sum(valid_pvalues) > 0) {
    results$p_adjusted <- NA_real_

    # Select correction method
    if (correction_method == "fdr") {
      # FDR control: BH or BY
      method <- if (hypotheses_independent) "BH" else "BY"
    } else {
      # FWER control: bonferroni, holm, hochberg, hommel
      method <- correction_method
    }

    results$p_adjusted[valid_pvalues] <- p.adjust(
      results$p_value[valid_pvalues],
      method = method
    )
    results$relevant <- results$p_adjusted < fdr_level
    results$correction_method <- method
  } else {
    results$p_adjusted <- NA_real_
    results$relevant <- FALSE
    results$correction_method <- NA_character_
  }

  # Constant features are never relevant
  results$relevant[results$type == "constant"] <- FALSE
  results$relevant[is.na(results$relevant)] <- FALSE

  # Sort by p-value
  results <- results[order(results$p_value, na.last = TRUE), ]
  rownames(results) <- NULL

  return(results)
}


#' Calculate Relevance for Multiclass Classification
#'
#' @inheritParams calculate_relevance_table
#' @param n_significant Minimum classes to be significant for
#'
#' @return data.frame with multiclass relevance information
#' @keywords internal
calculate_relevance_multiclass <- function(X, y,
                                          fdr_level = 0.05,
                                          hypotheses_independent = FALSE,
                                          correction_method = "fdr",
                                          n_significant = 1,
                                          test_binary_real = "mann") {

  classes <- unique(y)
  n_classes <- length(classes)

  if (n_classes <= 2) {
    warning("Only 2 or fewer classes found, using binary classification")
    return(calculate_relevance_table(
      X, y,
      fdr_level = fdr_level,
      hypotheses_independent = hypotheses_independent,
      correction_method = correction_method,
      ml_task = "classification",
      test_binary_real = test_binary_real
    ))
  }

  if (n_significant > n_classes) {
    stop("n_significant cannot exceed the number of classes")
  }

  # Test each class using one-vs-all approach
  class_results <- list()

  for (class_label in classes) {
    y_binary <- as.numeric(y == class_label)

    result <- calculate_relevance_table(
      X, y_binary,
      fdr_level = fdr_level,
      hypotheses_independent = hypotheses_independent,
      correction_method = correction_method,
      ml_task = "classification",
      test_binary_real = test_binary_real
    )

    class_results[[as.character(class_label)]] <- result
  }

  # Merge results from all classes
  relevance <- class_results[[1]][, c("feature", "type"), drop = FALSE]

  for (class_label in classes) {
    result <- class_results[[as.character(class_label)]]
    class_str <- as.character(class_label)

    relevance[[paste0("p_value_", class_str)]] <- result$p_value
    relevance[[paste0("p_adjusted_", class_str)]] <- result$p_adjusted
    relevance[[paste0("relevant_", class_str)]] <- result$relevant
  }

  # Count how many classes each feature is relevant for
  relevant_cols <- grep("^relevant_", names(relevance), value = TRUE)
  relevance$n_significant <- rowSums(relevance[, relevant_cols, drop = FALSE])
  relevance$relevant <- relevance$n_significant >= n_significant

  # Sort by number of significant classes (descending)
  relevance <- relevance[order(-relevance$n_significant), ]
  rownames(relevance) <- NULL

  return(relevance)
}


#' Classify Feature Type
#'
#' @param feature A numeric vector
#' @return "constant", "binary", or "real"
#' @keywords internal
classify_feature <- function(feature) {
  unique_vals <- unique(na.omit(feature))
  n_unique <- length(unique_vals)

  if (n_unique == 1) {
    return("constant")
  } else if (n_unique == 2) {
    return("binary")
  } else {
    return("real")
  }
}


#' Infer Machine Learning Task
#'
#' @param y Target vector
#' @return "classification" or "regression"
#' @keywords internal
infer_ml_task <- function(y) {
  if (is.factor(y) || is.logical(y) || is.character(y)) {
    return("classification")
  }

  # For numeric, check number of unique values
  n_unique <- length(unique(y))

  # If very few unique values (≤10) likely classification
  if (n_unique <= 10) {
    return("classification")
  }

  # If integer with moderate number of unique values (≤20) likely classification
  if (is.integer(y) && n_unique <= 20) {
    return("classification")
  }

  return("regression")
}


# ============================================================================
# Statistical Tests
# ============================================================================

#' Test Binary Target + Binary Feature
#'
#' Uses Fisher's exact test
#'
#' @param feature Binary feature vector
#' @param target Binary target vector
#' @return p-value
#' @keywords internal
test_binary_binary <- function(feature, target) {
  # Create contingency table
  tbl <- table(feature, target)

  # Fisher's exact test
  result <- fisher.test(tbl, alternative = "two.sided")

  return(result$p.value)
}


#' Test Binary Target + Real Feature
#'
#' Uses Mann-Whitney U test (default) or Kolmogorov-Smirnov test
#'
#' @param feature Real-valued feature vector
#' @param target Binary target vector (0/1)
#' @param test "mann" for Mann-Whitney U or "ks" for Kolmogorov-Smirnov
#' @return p-value
#' @keywords internal
test_binary_real <- function(feature, target, test = "mann") {
  # Split feature by target
  feature_0 <- feature[target == 0]
  feature_1 <- feature[target == 1]

  if (length(feature_0) == 0 || length(feature_1) == 0) {
    return(NA_real_)
  }

  if (test == "mann") {
    # Mann-Whitney U test (equivalent to Wilcoxon rank-sum)
    result <- wilcox.test(feature_1, feature_0, exact = FALSE, alternative = "two.sided")
  } else if (test == "ks") {
    # Kolmogorov-Smirnov test
    result <- ks.test(feature_1, feature_0, exact = FALSE, alternative = "two.sided")
  } else {
    stop("test must be 'mann' or 'ks'")
  }

  return(result$p.value)
}


#' Test Real Target + Binary Feature
#'
#' Uses Kolmogorov-Smirnov test
#'
#' @param feature Binary feature vector
#' @param target Real-valued target vector
#' @return p-value
#' @keywords internal
test_real_binary <- function(feature, target) {
  # Get unique feature values
  unique_vals <- unique(feature)

  if (length(unique_vals) != 2) {
    return(NA_real_)
  }

  # Split target by feature
  target_0 <- target[feature == unique_vals[1]]
  target_1 <- target[feature == unique_vals[2]]

  if (length(target_0) == 0 || length(target_1) == 0) {
    return(NA_real_)
  }

  # Kolmogorov-Smirnov test
  result <- ks.test(target_1, target_0, exact = FALSE, alternative = "two.sided")

  return(result$p.value)
}


#' Test Real Target + Real Feature
#'
#' Uses Kendall's tau correlation test
#'
#' @param feature Real-valued feature vector
#' @param target Real-valued target vector
#' @return p-value
#' @keywords internal
test_real_real <- function(feature, target) {
  # Kendall's tau correlation test
  result <- cor.test(feature, target, method = "kendall", exact = FALSE)

  return(result$p.value)
}


# ============================================================================
# Print and Summary Methods
# ============================================================================

#' @export
print.tsdiscov_feature_selection <- function(x, ...) {
  cat("Feature Selection Results\n")
  cat("==========================\n\n")
  cat(sprintf("Original features: %d\n", x$n_features_original))
  cat(sprintf("Selected features: %d (%.1f%%)\n",
              x$n_features_selected,
              100 * x$n_features_selected / x$n_features_original))
  cat(sprintf("FDR level: %.3f\n", x$fdr_level))
  cat(sprintf("ML task: %s\n\n", x$ml_task))

  if (x$n_features_selected > 0) {
    cat("Top 10 most relevant features:\n")
    top_features <- head(x$relevance_table[x$relevance_table$relevant, ], 10)
    print(top_features[, c("feature", "type", "p_value", "p_adjusted")], row.names = FALSE)
  } else {
    cat("No features passed the relevance threshold.\n")
    cat("Consider increasing fdr_level or checking your data.\n")
  }

  invisible(x)
}


#' @export
summary.tsdiscov_feature_selection <- function(object, ...) {
  cat("Feature Selection Summary\n")
  cat("=========================\n\n")

  cat(sprintf("Total features: %d\n", object$n_features_original))
  cat(sprintf("Relevant features: %d (%.1f%%)\n",
              object$n_features_selected,
              100 * object$n_features_selected / object$n_features_original))
  cat(sprintf("FDR level: %.3f\n", object$fdr_level))
  cat(sprintf("ML task: %s\n\n", object$ml_task))

  # Feature type breakdown
  type_counts <- table(object$relevance_table$type)
  cat("Features by type:\n")
  print(type_counts)
  cat("\n")

  # Relevant features by type
  relevant_table <- object$relevance_table[object$relevance_table$relevant, ]
  if (nrow(relevant_table) > 0) {
    relevant_type_counts <- table(relevant_table$type)
    cat("Relevant features by type:\n")
    print(relevant_type_counts)
    cat("\n")
  }

  # P-value distribution
  valid_pvals <- object$relevance_table$p_value[!is.na(object$relevance_table$p_value)]
  if (length(valid_pvals) > 0) {
    cat("P-value distribution:\n")
    cat(sprintf("  Min: %.4f\n", min(valid_pvals)))
    cat(sprintf("  Q1:  %.4f\n", quantile(valid_pvals, 0.25)))
    cat(sprintf("  Median: %.4f\n", median(valid_pvals)))
    cat(sprintf("  Q3:  %.4f\n", quantile(valid_pvals, 0.75)))
    cat(sprintf("  Max: %.4f\n", max(valid_pvals)))
  }

  invisible(object)
}
