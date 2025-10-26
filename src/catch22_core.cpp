#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <map>
using namespace Rcpp;

// Helper function: Calculate histogram mode
double histogram_mode(NumericVector x, int num_bins) {
    if (x.size() == 0) return NA_REAL;

    double min_val = *std::min_element(x.begin(), x.end());
    double max_val = *std::max_element(x.begin(), x.end());

    if (max_val == min_val) return min_val;

    std::vector<int> bins(num_bins, 0);
    double bin_width = (max_val - min_val) / num_bins;

    // Fill histogram
    for (int i = 0; i < x.size(); i++) {
        int bin = static_cast<int>((x[i] - min_val) / bin_width);
        if (bin >= num_bins) bin = num_bins - 1;
        bins[bin]++;
    }

    // Find mode
    int max_count = *std::max_element(bins.begin(), bins.end());
    double mode_count = static_cast<double>(max_count) / x.size();

    return mode_count;
}

// [[Rcpp::export]]
double cpp_histogram_mode_5(NumericVector x) {
    return histogram_mode(x, 5);
}

// [[Rcpp::export]]
double cpp_histogram_mode_10(NumericVector x) {
    return histogram_mode(x, 10);
}

// [[Rcpp::export]]
double cpp_outlier_timing(NumericVector x, bool positive = true, double threshold = 0.01) {
    int n = x.size();
    if (n < 3) return NA_REAL;

    // Calculate median absolute deviation
    NumericVector sorted = clone(x);
    std::sort(sorted.begin(), sorted.end());
    double median = sorted[n / 2];

    NumericVector abs_dev(n);
    for (int i = 0; i < n; i++) {
        abs_dev[i] = std::abs(x[i] - median);
    }
    std::sort(abs_dev.begin(), abs_dev.end());
    double mad = abs_dev[n / 2];

    if (mad == 0) return NA_REAL;

    // Find outliers
    std::vector<int> outlier_positions;
    double z_threshold = 3.0; // Standard threshold for outliers

    for (int i = 0; i < n; i++) {
        double z_score = std::abs(x[i] - median) / (1.4826 * mad);

        if (positive && x[i] > median && z_score > z_threshold) {
            outlier_positions.push_back(i);
        } else if (!positive && x[i] < median && z_score > z_threshold) {
            outlier_positions.push_back(i);
        }
    }

    if (outlier_positions.empty()) return 0.0;

    // Calculate mean relative position
    double mean_pos = 0.0;
    for (size_t i = 0; i < outlier_positions.size(); i++) {
        mean_pos += static_cast<double>(outlier_positions[i]) / n;
    }
    return mean_pos / outlier_positions.size();
}

// [[Rcpp::export]]
double cpp_outlier_timing_pos(NumericVector x) {
    return cpp_outlier_timing(x, true, 0.01);
}

// [[Rcpp::export]]
double cpp_outlier_timing_neg(NumericVector x) {
    return cpp_outlier_timing(x, false, 0.01);
}

// [[Rcpp::export]]
double cpp_long_stretch(NumericVector x, bool above_mean = true) {
    int n = x.size();
    if (n == 0) return NA_REAL;

    // Calculate mean
    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= n;

    // Find longest stretch
    int max_stretch = 0;
    int current_stretch = 0;

    for (int i = 0; i < n; i++) {
        if ((above_mean && x[i] > mean) || (!above_mean && x[i] < mean)) {
            current_stretch++;
            if (current_stretch > max_stretch) {
                max_stretch = current_stretch;
            }
        } else {
            current_stretch = 0;
        }
    }

    return static_cast<double>(max_stretch);
}

// [[Rcpp::export]]
double cpp_stretch_high(NumericVector x) {
    return cpp_long_stretch(x, true);
}

// [[Rcpp::export]]
double cpp_stretch_low(NumericVector x) {
    return cpp_long_stretch(x, false);
}

// [[Rcpp::export]]
double cpp_transition_variance(NumericVector x, int num_bins = 3) {
    int n = x.size();
    if (n < 2) return NA_REAL;

    // Discretize into bins
    double min_val = *std::min_element(x.begin(), x.end());
    double max_val = *std::max_element(x.begin(), x.end());

    if (max_val == min_val) return 0.0;

    std::vector<int> discrete(n);
    double bin_width = (max_val - min_val) / num_bins;

    for (int i = 0; i < n; i++) {
        discrete[i] = static_cast<int>((x[i] - min_val) / bin_width);
        if (discrete[i] >= num_bins) discrete[i] = num_bins - 1;
    }

    // Build transition matrix
    std::vector<std::vector<int>> trans_matrix(num_bins, std::vector<int>(num_bins, 0));
    for (int i = 0; i < n - 1; i++) {
        trans_matrix[discrete[i]][discrete[i + 1]]++;
    }

    // Calculate variance of transition matrix diagonal
    std::vector<double> diag(num_bins);
    int total = 0;
    for (int i = 0; i < num_bins; i++) {
        for (int j = 0; j < num_bins; j++) {
            total += trans_matrix[i][j];
        }
    }

    if (total == 0) return 0.0;

    // Normalize and get diagonal
    for (int i = 0; i < num_bins; i++) {
        int row_sum = 0;
        for (int j = 0; j < num_bins; j++) {
            row_sum += trans_matrix[i][j];
        }
        diag[i] = row_sum > 0 ? static_cast<double>(trans_matrix[i][i]) / row_sum : 0.0;
    }

    // Calculate variance
    double mean = 0.0;
    for (int i = 0; i < num_bins; i++) mean += diag[i];
    mean /= num_bins;

    double var = 0.0;
    for (int i = 0; i < num_bins; i++) {
        var += (diag[i] - mean) * (diag[i] - mean);
    }
    return var / num_bins;
}
