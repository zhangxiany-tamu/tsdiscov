#' Validate and convert multivariate input to matrix format
#'
#' @param X Either a matrix (N x T) or list of N numeric vectors
#' @return Matrix with N rows (series) and T columns (timepoints)
#' @keywords internal
validate_multivariate_input <- function(X) {
  # Case 1: Already a matrix
  if (is.matrix(X)) {
    if (!is.numeric(X)) {
      stop("Matrix X must be numeric")
    }
    if (nrow(X) < 2) {
      stop("Need at least 2 time series (N >= 2)")
    }
    if (ncol(X) < 3) {
      stop("Each time series must have at least 3 observations (T >= 3)")
    }
    return(X)
  }

  # Case 2: List of vectors
  if (is.list(X)) {
    # Check all elements are numeric vectors
    if (!all(sapply(X, is.numeric))) {
      stop("All elements in list X must be numeric vectors")
    }

    # Check lengths
    lengths <- sapply(X, length)
    if (length(unique(lengths)) > 1) {
      stop("All time series must have equal length")
    }

    N <- length(X)
    T_len <- lengths[1]

    if (N < 2) {
      stop("Need at least 2 time series (N >= 2)")
    }
    if (T_len < 3) {
      stop("Each time series must have at least 3 observations (T >= 3)")
    }

    # Convert to matrix (each row is a series)
    X_matrix <- do.call(rbind, X)
    return(X_matrix)
  }

  # Case 3: Data frame - convert to matrix
  if (is.data.frame(X)) {
    if (ncol(X) < 2) {
      stop("Data frame must have at least 2 columns (series)")
    }
    X_matrix <- as.matrix(X)
    # Transpose so rows are series
    return(t(X_matrix))
  }

  stop("X must be a matrix (N x T), list of N vectors, or data.frame")
}

#' Standardize multivariate time series
#'
#' Standardizes each series (row) to have mean 0 and sd 1
#'
#' @param X Matrix (N x T)
#' @return Standardized matrix
#' @keywords internal
standardize_multivariate <- function(X) {
  # Standardize each row (series)
  X_std <- t(scale(t(X)))

  # Handle constant series (sd = 0)
  constant_rows <- apply(X, 1, function(row) sd(row) == 0)
  if (any(constant_rows)) {
    warning(sprintf("%d series are constant and will be set to zero after standardization",
                    sum(constant_rows)))
    X_std[constant_rows, ] <- 0
  }

  return(X_std)
}

#' Get names of multivariate feature sets
#'
#' @return Character vector of available feature set names
#' @keywords internal
get_multivariate_feature_sets <- function() {
  c("pca", "correlation", "covariance", "sync", "diversity")
}
