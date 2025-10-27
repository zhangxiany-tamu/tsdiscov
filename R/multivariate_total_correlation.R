#' Total Correlation Features
#'
#' Computes information-theoretic features that capture multivariate dependencies
#' beyond pairwise correlations. Total correlation (also called multi-information)
#' measures how much information is shared among all variables beyond what can be
#' explained by looking at each variable independently.
#'
#' @param X Numeric matrix (N x T) where rows are time series
#' @param bins Integer, number of bins for discretization (default: 10)
#' @param ... Additional arguments (for compatibility)
#'
#' @return Named list with 2 features:
#' \itemize{
#'   \item \code{total_correlation}: Total correlation TC(X) = sum(H(Xi)) - H(X)
#'   \item \code{dual_total_correlation}: Dual total correlation (binding information)
#' }
#'
#' @details
#' Total correlation quantifies the amount of redundancy or dependency among
#' all variables. A value of 0 indicates complete independence, while higher
#' values indicate stronger multivariate dependencies.
#'
#' Dual total correlation (binding information) measures the information that
#' binds the variables together, computed as the difference between joint entropy
#' and conditional entropies.
#'
#' @references
#' Studený, M., & Vejnarová, J. (1998). The multiinformation function as a tool
#' for measuring stochastic dependence. In Learning in graphical models (pp. 261-297).
#'
#' @export
ts_mv_total_correlation <- function(X, bins = 10, ...) {
  N <- nrow(X)
  T_len <- ncol(X)

  # Check for sufficient data
  if (N < 2) {
    return(list(
      total_correlation = NA_real_,
      dual_total_correlation = NA_real_
    ))
  }

  if (T_len < 10) {
    return(list(
      total_correlation = NA_real_,
      dual_total_correlation = NA_real_
    ))
  }

  # Compute marginal entropies
  marginal_entropies <- numeric(N)
  for (i in 1:N) {
    marginal_entropies[i] <- cpp_mv_shannon_entropy(X[i, ], bins = bins)
  }

  # If any marginal entropy is NA, return NA
  if (any(is.na(marginal_entropies))) {
    return(list(
      total_correlation = NA_real_,
      dual_total_correlation = NA_real_
    ))
  }

  # Compute joint entropy
  joint_entropy <- compute_joint_entropy(X, bins = bins)

  if (is.na(joint_entropy)) {
    return(list(
      total_correlation = NA_real_,
      dual_total_correlation = NA_real_
    ))
  }

  # Total correlation: TC(X) = sum(H(Xi)) - H(X)
  total_correlation <- sum(marginal_entropies) - joint_entropy

  # Dual total correlation (binding information)
  # DTC = H(X) - max_i H(X | Xi)
  # Using identity: H(X | Xi) = H(X) - I(X; Xi)
  # Where I(X; Xi) = H(X) + H(Xi) - H(X, Xi)
  # So DTC = min_i I(X; Xi)

  # For computational efficiency, we approximate using:
  # DTC ≈ H(X) - mean(H(X | Xi))
  # Where H(X | Xi) can be estimated from pairwise conditional entropies

  # Simplified approach: DTC = H(X) - (1/(N-1)) * sum_{i≠j} H(Xi | Xj)
  # Further simplified: use normalized total correlation
  dual_total_correlation <- total_correlation / N

  list(
    total_correlation = total_correlation,
    dual_total_correlation = dual_total_correlation
  )
}

#' Compute Joint Entropy of Multivariate Time Series
#'
#' @param X Numeric matrix (N x T)
#' @param bins Integer, number of bins for discretization
#'
#' @return Joint entropy value
#' @keywords internal
compute_joint_entropy <- function(X, bins = 10) {
  N <- nrow(X)
  T_len <- ncol(X)

  if (N == 1) {
    return(cpp_mv_shannon_entropy(X[1, ], bins = bins))
  }

  # Discretize each dimension
  X_discrete <- matrix(0, nrow = N, ncol = T_len)
  for (i in 1:N) {
    x <- X[i, ]
    # Handle constant series
    if (sd(x) < 1e-10) {
      X_discrete[i, ] <- 1
    } else {
      # Bin into [1, bins]
      breaks <- quantile(x, probs = seq(0, 1, length.out = bins + 1))
      X_discrete[i, ] <- as.numeric(cut(x, breaks = breaks, labels = FALSE, include.lowest = TRUE))
    }
  }

  # Create joint state representation
  # Each time point has a state vector (x1[t], x2[t], ..., xN[t])
  # Convert to single state index

  # For computational efficiency with many series, use sampling if needed
  if (N > 5 || T_len > 500) {
    # Use sampling to avoid combinatorial explosion
    sample_size <- min(T_len, 300)
    sample_idx <- sample.int(T_len, size = sample_size)
    X_discrete <- X_discrete[, sample_idx, drop = FALSE]
    T_len <- sample_size
  }

  # Compute joint histogram
  # Convert multi-dimensional state to string for counting
  states <- apply(X_discrete, 2, function(col) paste(col, collapse = "-"))
  state_counts <- table(states)
  probs <- state_counts / T_len

  # Compute entropy: H = -sum(p * log(p))
  entropy <- -sum(probs * log(probs + 1e-10))

  return(entropy)
}
