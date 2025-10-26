#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <map>
using namespace Rcpp;

// Recurrence Quantification Analysis (RQA)
// Analyzes recurrence plot to extract nonlinear dynamics features
// [[Rcpp::export]]
List cpp_recurrence_analysis(NumericVector x, double threshold_percent = 0.1, int min_line_length = 2, int max_points = 2000) {
    int n = x.size();

    // Default result
    List default_result = List::create(
        Named("recurrence_rate") = NA_REAL,
        Named("determinism") = NA_REAL,
        Named("laminarity") = NA_REAL,
        Named("longest_diagonal") = NA_REAL,
        Named("entropy_diagonal") = NA_REAL
    );

    // Need sufficient data
    if (n < 50) return default_result;

    // Sample if series is too long (for computational efficiency)
    std::vector<double> y;
    int sample_rate = 1;
    if (n > max_points) {
        sample_rate = n / max_points;
        for (int i = 0; i < n; i += sample_rate) {
            y.push_back(x[i]);
        }
    } else {
        for (int i = 0; i < n; i++) {
            y.push_back(x[i]);
        }
    }

    int m = y.size();
    if (m < 50) return default_result;

    // Calculate threshold based on percentile of distances
    std::vector<double> distances;
    distances.reserve(m * 10);  // Sample some distances

    for (int i = 0; i < std::min(m, 200); i += 5) {
        for (int j = i + 1; j < std::min(m, 200); j += 5) {
            double dist = std::abs(y[i] - y[j]);
            distances.push_back(dist);
        }
    }

    if (distances.empty()) return default_result;

    std::sort(distances.begin(), distances.end());
    int threshold_idx = static_cast<int>(distances.size() * threshold_percent);
    double threshold = distances[threshold_idx];

    // Build recurrence matrix (sparse representation)
    // Only store recurrent points
    std::vector<std::vector<int>> recurrence(m);

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            if (std::abs(y[i] - y[j]) < threshold) {
                recurrence[i].push_back(j);
            }
        }
    }

    // Calculate Recurrence Rate (RR)
    int total_recurrent = 0;
    for (int i = 0; i < m; i++) {
        total_recurrent += recurrence[i].size();
    }
    double recurrence_rate = static_cast<double>(total_recurrent) / (m * m);

    // Find diagonal lines (determinism indicator)
    std::vector<int> diagonal_lengths;

    for (int k = 1 - m; k < m; k++) {
        int current_length = 0;

        for (int i = 0; i < m; i++) {
            int j = i + k;
            if (j >= 0 && j < m) {
                // Check if (i, j) is recurrent
                bool is_recurrent = false;
                for (int rec_j : recurrence[i]) {
                    if (rec_j == j) {
                        is_recurrent = true;
                        break;
                    }
                }

                if (is_recurrent) {
                    current_length++;
                } else {
                    if (current_length >= min_line_length) {
                        diagonal_lengths.push_back(current_length);
                    }
                    current_length = 0;
                }
            }
        }

        if (current_length >= min_line_length) {
            diagonal_lengths.push_back(current_length);
        }
    }

    // Calculate Determinism (DET)
    int diagonal_points = 0;
    for (int len : diagonal_lengths) {
        diagonal_points += len;
    }

    double determinism = (total_recurrent > 0) ?
        static_cast<double>(diagonal_points) / total_recurrent : 0.0;

    // Find vertical lines (laminarity indicator)
    std::vector<int> vertical_lengths;

    for (int j = 0; j < m; j++) {
        int current_length = 0;

        for (int i = 0; i < m; i++) {
            // Check if (i, j) is recurrent
            bool is_recurrent = false;
            for (int rec_j : recurrence[i]) {
                if (rec_j == j) {
                    is_recurrent = true;
                    break;
                }
            }

            if (is_recurrent) {
                current_length++;
            } else {
                if (current_length >= min_line_length) {
                    vertical_lengths.push_back(current_length);
                }
                current_length = 0;
            }
        }

        if (current_length >= min_line_length) {
            vertical_lengths.push_back(current_length);
        }
    }

    // Calculate Laminarity (LAM)
    int vertical_points = 0;
    for (int len : vertical_lengths) {
        vertical_points += len;
    }

    double laminarity = (total_recurrent > 0) ?
        static_cast<double>(vertical_points) / total_recurrent : 0.0;

    // Longest diagonal line
    int longest_diagonal = 0;
    if (!diagonal_lengths.empty()) {
        longest_diagonal = *std::max_element(diagonal_lengths.begin(), diagonal_lengths.end());
    }

    // Shannon entropy of diagonal line length distribution
    double entropy_diagonal = 0.0;
    if (!diagonal_lengths.empty()) {
        std::map<int, int> length_counts;
        for (int len : diagonal_lengths) {
            length_counts[len]++;
        }

        int total_lines = diagonal_lengths.size();
        for (const auto& pair : length_counts) {
            double p = static_cast<double>(pair.second) / total_lines;
            entropy_diagonal -= p * std::log(p);
        }
    }

    return List::create(
        Named("recurrence_rate") = recurrence_rate,
        Named("determinism") = determinism,
        Named("laminarity") = laminarity,
        Named("longest_diagonal") = static_cast<double>(longest_diagonal),
        Named("entropy_diagonal") = entropy_diagonal
    );
}
