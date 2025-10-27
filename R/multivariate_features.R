#' Extract features from multivariate time series
#'
#' Extracts features that characterize the joint behavior of multiple time series.
#' Supports analysis of sensor arrays, multi-channel signals, and other multivariate
#' engineering data.
#'
#' @param X Either a matrix (N x T) where N is number of series and T is timepoints,
#'          or a list of N numeric vectors of equal length T
#' @param features Character vector of feature sets to extract. Options:
#'   \itemize{
#'     \item "all" - All available multivariate features (default)
#'     \item "pca" - PCA-based features (15 features)
#'     \item "correlation" - Correlation structure features (16 features)
#'     \item "covariance" - Covariance matrix features (7 features)
#'     \item "spectral" - Cross-spectral and coherence features (10 features)
#'     \item "sync" - Synchronization features (10 features)
#'     \item "diversity" - Diversity and complexity features (8 features)
#'   }
#'   Can specify multiple sets, e.g., c("pca", "correlation")
#' @param standardize Logical; if TRUE (default), standardizes each series to
#'   have mean 0 and standard deviation 1 before computing features
#' @return Named list of multivariate features. The number of features depends
#'   on which feature sets are requested.
#' @export
#' @examples
#' \dontrun{
#' # Example 1: Matrix input (5 sensors, 200 timepoints)
#' set.seed(123)
#' X_matrix <- matrix(rnorm(5 * 200), nrow = 5, ncol = 200)
#' # Add correlation structure
#' X_matrix[2,] <- 0.7 * X_matrix[1,] + 0.3 * rnorm(200)
#' features <- ts_features_multivariate(X_matrix)
#'
#' # Example 2: List input (3 channels)
#' X_list <- list(
#'   channel1 = rnorm(150),
#'   channel2 = rnorm(150),
#'   channel3 = rnorm(150)
#' )
#' features <- ts_features_multivariate(X_list, standardize = TRUE)
#'
#' # Example 3: Extract specific feature sets
#' pca_feats <- ts_features_multivariate(X_matrix, features = "pca")
#' corr_feats <- ts_features_multivariate(X_matrix, features = c("pca", "correlation"))
#' }
ts_features_multivariate <- function(X,
                                     features = "all",
                                     standardize = TRUE) {
  # Validate and convert input to matrix format
  X <- validate_multivariate_input(X)

  # Check for missing data
  if (anyNA(X)) {
    stop("Missing data not supported. Please provide complete data.")
  }

  # Standardize if requested
  if (standardize) {
    X <- standardize_multivariate(X)
  }

  # Expand "all" to all available feature sets
  if (length(features) == 1 && features == "all") {
    features <- get_multivariate_feature_sets()
  }

  # Validate feature set names
  valid_sets <- get_multivariate_feature_sets()
  invalid <- setdiff(features, valid_sets)
  if (length(invalid) > 0) {
    stop(sprintf("Invalid feature set(s): %s\nValid options: %s",
                 paste(invalid, collapse = ", "),
                 paste(valid_sets, collapse = ", ")))
  }

  # Extract requested features
  result <- list()

  if ("pca" %in% features) {
    result <- c(result, ts_mv_pca(X))
  }

  if ("correlation" %in% features) {
    result <- c(result, ts_mv_correlation(X))
  }

  if ("covariance" %in% features) {
    result <- c(result, ts_mv_covariance(X))
  }

  if ("sync" %in% features) {
    result <- c(result, ts_mv_synchronization(X))
  }

  if ("diversity" %in% features) {
    result <- c(result, ts_mv_diversity(X))
  }

  return(result)
}

#' Extract multivariate features as data frame
#'
#' Wrapper around ts_features_multivariate() that returns results as a 1-row
#' data frame. This format is convenient for batch processing of multiple
#' multivariate datasets.
#'
#' @param X Either a matrix (N x T) or list of N numeric vectors
#' @param features Character vector of feature sets to extract (default "all")
#' @param standardize Logical; standardize series before extraction (default TRUE)
#' @return A 1-row data.frame with one column per feature
#' @export
#' @examples
#' \dontrun{
#' # Extract features from multiple multivariate datasets
#' datasets <- list(
#'   system1 = matrix(rnorm(3 * 100), nrow = 3),
#'   system2 = matrix(rnorm(3 * 100), nrow = 3),
#'   system3 = matrix(rnorm(3 * 100), nrow = 3)
#' )
#'
#' results <- do.call(rbind, lapply(datasets, ts_features_multivariate_df))
#' rownames(results) <- names(datasets)
#' }
ts_features_multivariate_df <- function(X,
                                        features = "all",
                                        standardize = TRUE) {
  result <- ts_features_multivariate(X, features = features, standardize = standardize)
  as.data.frame(t(unlist(result)))
}
