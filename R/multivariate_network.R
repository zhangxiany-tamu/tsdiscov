#' Graph/Network Topology Features
#'
#' Constructs a correlation network where nodes represent time series and edges
#' represent strong correlations, then computes standard network topology metrics.
#' These features capture connectivity patterns, community structure, and the
#' overall network architecture of the multivariate system.
#'
#' @param X Numeric matrix (N x T) where rows are time series
#' @param threshold Numeric; correlation threshold for edge creation (default: 0.5).
#'   Only correlations with |r| >= threshold create edges.
#' @param ... Additional arguments (for compatibility)
#'
#' @return Named list with 5 features:
#' \itemize{
#'   \item \code{network_density}: Fraction of possible edges that exist
#'   \item \code{network_clustering}: Average clustering coefficient
#'   \item \code{network_assortativity}: Degree assortativity (degree correlation)
#'   \item \code{network_modularity}: Newman modularity (community structure strength)
#'   \item \code{network_centralization}: Network centralization (hub dominance)
#' }
#'
#' @details
#' The correlation network is constructed by:
#' \enumerate{
#'   \item Computing pairwise correlations between all series
#'   \item Creating edges for correlations with |r| >= threshold
#'   \item Computing network metrics on the resulting graph
#' }
#'
#' Network metrics interpretation:
#' \itemize{
#'   \item **Density**: Ranges from 0 (no edges) to 1 (complete graph). High density
#'         indicates many strong pairwise relationships.
#'   \item **Clustering**: Measures local connectivity. High clustering indicates
#'         tightly connected groups.
#'   \item **Assortativity**: Measures if high-degree nodes connect to other
#'         high-degree nodes. Positive = assortative, negative = disassortative.
#'   \item **Modularity**: Measures community structure strength. Higher values
#'         indicate distinct modules/communities.
#'   \item **Centralization**: Measures if network has dominant hub nodes. High
#'         values indicate star-like topology.
#' }
#'
#' @references
#' Newman, M. E. (2003). The structure and function of complex networks.
#' SIAM review, 45(2), 167-256.
#'
#' @export
ts_mv_network <- function(X, threshold = 0.5, ...) {
  N <- nrow(X)

  # Check for sufficient data
  if (N < 3) {
    return(list(
      network_density = NA_real_,
      network_clustering = NA_real_,
      network_assortativity = NA_real_,
      network_modularity = NA_real_,
      network_centralization = NA_real_
    ))
  }

  # Compute correlation matrix
  R <- cor(t(X))

  # Create adjacency matrix (undirected graph)
  # Edge exists if |correlation| >= threshold
  A <- abs(R) >= threshold
  diag(A) <- FALSE  # Remove self-loops

  # Convert to numeric for calculations
  A_num <- A * 1.0

  # Compute network metrics

  # 1. Density: fraction of existing edges
  n_possible_edges <- N * (N - 1) / 2
  n_edges <- sum(A_num) / 2  # Divide by 2 for undirected
  network_density <- n_edges / n_possible_edges

  # 2. Clustering coefficient: average over all nodes
  clustering_coeffs <- numeric(N)
  for (i in 1:N) {
    neighbors <- which(A[i, ])
    k_i <- length(neighbors)

    if (k_i < 2) {
      clustering_coeffs[i] <- 0
    } else {
      # Count edges among neighbors
      subgraph <- A[neighbors, neighbors]
      n_neighbor_edges <- sum(subgraph) / 2
      n_possible_neighbor_edges <- k_i * (k_i - 1) / 2
      clustering_coeffs[i] <- n_neighbor_edges / n_possible_neighbor_edges
    }
  }
  network_clustering <- mean(clustering_coeffs)

  # 3. Degree assortativity: correlation of degrees at edge endpoints
  degrees <- rowSums(A_num)

  if (n_edges == 0) {
    network_assortativity <- 0
  } else {
    # Get degree pairs for all edges
    edge_list <- which(A & upper.tri(A), arr.ind = TRUE)
    if (nrow(edge_list) > 0) {
      degree_i <- degrees[edge_list[, 1]]
      degree_j <- degrees[edge_list[, 2]]

      # Check for zero variance or NA (all degrees same or too few edges)
      sd_i <- sd(degree_i)
      sd_j <- sd(degree_j)

      if (is.na(sd_i) || is.na(sd_j) || sd_i < 1e-10 || sd_j < 1e-10) {
        network_assortativity <- 0
      } else {
        network_assortativity <- cor(degree_i, degree_j)
      }
    } else {
      network_assortativity <- 0
    }
  }

  # 4. Modularity: use simple greedy community detection
  # For computational efficiency, use a simplified modularity calculation
  # based on degree-based null model
  if (n_edges == 0) {
    network_modularity <- 0
  } else {
    # Simple community detection: group by degree similarity
    # This is a fast approximation
    degree_order <- order(degrees, decreasing = TRUE)
    n_communities <- min(max(2, floor(N / 3)), N)
    community <- rep(1:n_communities, length.out = N)[rank(-degrees, ties.method = "first")]

    # Compute modularity
    m <- n_edges
    Q <- 0
    for (c in unique(community)) {
      nodes_in_c <- which(community == c)
      for (i in nodes_in_c) {
        for (j in nodes_in_c) {
          if (i < j) {
            A_ij <- A_num[i, j]
            k_i <- degrees[i]
            k_j <- degrees[j]
            expected <- (k_i * k_j) / (2 * m)
            Q <- Q + (A_ij - expected)
          }
        }
      }
    }
    network_modularity <- Q / (2 * m)
  }

  # 5. Centralization: based on degree centrality
  # Star centralization formula
  max_degree <- max(degrees)
  if (max_degree == 0) {
    network_centralization <- 0
  } else {
    sum_diffs <- sum(max_degree - degrees)
    max_sum_diffs <- (N - 1) * (N - 2)  # Max for star graph
    network_centralization <- if (max_sum_diffs > 0) sum_diffs / max_sum_diffs else 0
  }

  list(
    network_density = network_density,
    network_clustering = network_clustering,
    network_assortativity = ifelse(is.na(network_assortativity), 0, network_assortativity),
    network_modularity = network_modularity,
    network_centralization = network_centralization
  )
}
