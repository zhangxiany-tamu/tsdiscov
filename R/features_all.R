#' Extract all time series features (univariate and/or multivariate)
#'
#' Unified interface for extracting comprehensive time series features.
#' Automatically detects input type and extracts appropriate features.
#'
#' @param X Numeric vector (univariate), matrix (N x T), or list of N vectors
#' @param feature_type Character string specifying which features to extract:
#'   \itemize{
#'     \item "auto" - Automatically detect (default): univariate for vectors,
#'           multivariate for matrices/lists with N >= 2
#'     \item "univariate" - Extract only univariate features
#'     \item "multivariate" - Extract only multivariate features (requires N >= 2)
#'     \item "both" - Extract both univariate and multivariate features
#'   }
#' @param univariate_summary Character string for multivariate input:
#'   \itemize{
#'     \item "none" - No univariate features (default for multivariate)
#'     \item "aggregate" - Aggregate univariate features across series
#'           (mean, median, max, min, sd of each feature)
#'     \item "per_series" - Individual univariate features for each series
#'           (can produce many features for large N)
#'   }
#' @param multivariate_sets Character vector of multivariate feature sets
#'   (default "all"). Options: "all", "pca", "correlation", "covariance",
#'   "sync", "spectral", "diversity", or combinations
#' @param standardize Logical; standardize multivariate series (default TRUE)
#' @param ... Additional arguments passed to feature extraction functions
#' @return Named list of features. Length depends on feature_type and input.
#' @export
#' @examples
#' \dontrun{
#' # Univariate: single time series
#' x <- rnorm(200)
#' features <- ts_features_all(x)
#' length(features)  # 352 univariate features
#'
#' # Multivariate: multiple time series
#' X <- matrix(rnorm(5 * 200), nrow = 5, ncol = 200)
#' features <- ts_features_all(X)
#' length(features)  # 66 multivariate features
#'
#' # Both univariate and multivariate
#' features <- ts_features_all(X, feature_type = "both",
#'                              univariate_summary = "aggregate")
#' length(features)  # 66 + aggregated univariate
#'
#' # Only multivariate sync and spectral features
#' features <- ts_features_all(X, feature_type = "multivariate",
#'                              multivariate_sets = c("sync", "spectral"))
#' }
ts_features_all <- function(X,
                            feature_type = c("auto", "univariate", "multivariate", "both"),
                            univariate_summary = c("none", "aggregate", "per_series"),
                            multivariate_sets = "all",
                            standardize = TRUE,
                            ...) {
  feature_type <- match.arg(feature_type)
  univariate_summary <- match.arg(univariate_summary)

  # Determine if input is univariate or multivariate
  is_univariate <- (is.vector(X) && !is.list(X)) || (is.matrix(X) && nrow(X) == 1)

  # Auto-detect feature type
  if (feature_type == "auto") {
    feature_type <- if (is_univariate) "univariate" else "multivariate"
  }

  # Validate input
  if (feature_type == "multivariate" && is_univariate) {
    stop("Cannot extract multivariate features from single time series. ",
         "Use feature_type = 'univariate' or provide multiple series.")
  }

  if (feature_type == "both" && is_univariate) {
    warning("Single series provided. Returning univariate features only.")
    feature_type <- "univariate"
  }

  # Extract features based on type
  result <- list()

  # ========================================================================
  # Univariate features
  # ========================================================================
  if (feature_type %in% c("univariate", "both")) {

    if (is_univariate) {
      # Single series: extract univariate features
      result <- ts_features(as.numeric(X), ...)

    } else {
      # Multiple series: handle based on summary type
      if (univariate_summary == "none") {
        # Don't extract univariate features
        # (This only makes sense if feature_type = "both")

      } else if (univariate_summary == "per_series") {
        # Extract features for each series individually
        N <- if (is.matrix(X)) nrow(X) else length(X)

        for (i in 1:N) {
          series <- if (is.matrix(X)) X[i,] else X[[i]]
          series_features <- ts_features(as.numeric(series), ...)

          # Prefix with series number
          names(series_features) <- paste0("series", i, "_", names(series_features))
          result <- c(result, series_features)
        }

      } else if (univariate_summary == "aggregate") {
        # Extract features for each series, then aggregate
        result <- c(result, aggregate_univariate_features(X, ...))
      }
    }
  }

  # ========================================================================
  # Multivariate features
  # ========================================================================
  if (feature_type %in% c("multivariate", "both")) {

    if (!is_univariate) {
      mv_features <- ts_features_multivariate(
        X,
        features = multivariate_sets,
        standardize = standardize
      )
      result <- c(result, mv_features)
    }
  }

  return(result)
}

#' Aggregate univariate features across multiple series
#'
#' Computes summary statistics (mean, median, max, min, sd) of univariate
#' features across all series in a multivariate dataset.
#'
#' @param X Matrix (N x T) or list of N numeric vectors
#' @param ... Additional arguments passed to ts_features()
#' @return Named list of aggregated features
#' @keywords internal
aggregate_univariate_features <- function(X, ...) {
  # Convert to matrix if needed
  if (is.list(X) && !is.data.frame(X)) {
    X <- do.call(rbind, X)
  }

  N <- nrow(X)

  # Extract features for each series
  feature_matrix <- matrix(NA, nrow = N, ncol = 0)
  feature_names <- NULL

  for (i in 1:N) {
    series_features <- ts_features(as.numeric(X[i,]), ...)

    if (i == 1) {
      feature_names <- names(series_features)
      feature_matrix <- matrix(NA, nrow = N, ncol = length(feature_names))
      colnames(feature_matrix) <- feature_names
    }

    feature_matrix[i, ] <- unlist(series_features)
  }

  # Compute aggregated statistics
  result <- list()

  for (feat_name in feature_names) {
    feat_values <- feature_matrix[, feat_name]

    # Only compute if we have at least one non-NA value
    if (sum(!is.na(feat_values)) > 0) {
      result[[paste0(feat_name, "_mean")]] <- mean(feat_values, na.rm = TRUE)
      result[[paste0(feat_name, "_max")]] <- suppressWarnings(max(feat_values, na.rm = TRUE))
      result[[paste0(feat_name, "_min")]] <- suppressWarnings(min(feat_values, na.rm = TRUE))
      result[[paste0(feat_name, "_sd")]] <- sd(feat_values, na.rm = TRUE)
    } else {
      # All NA - return NA for all aggregations
      result[[paste0(feat_name, "_mean")]] <- NA_real_
      result[[paste0(feat_name, "_max")]] <- NA_real_
      result[[paste0(feat_name, "_min")]] <- NA_real_
      result[[paste0(feat_name, "_sd")]] <- NA_real_
    }
  }

  return(result)
}

#' Extract all features as data frame
#'
#' Wrapper around ts_features_all() that returns a 1-row data frame.
#' Convenient for batch processing.
#'
#' @param X Numeric vector, matrix, or list
#' @param ... Arguments passed to ts_features_all()
#' @return 1-row data.frame with features as columns
#' @export
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(5 * 200), nrow = 5)
#' df <- ts_features_all_df(X)
#' dim(df)  # 1 x 66
#' }
ts_features_all_df <- function(X, ...) {
  result <- ts_features_all(X, ...)
  as.data.frame(t(unlist(result)))
}
