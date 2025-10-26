#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
double cpp_automutual_info_lag2(NumericVector x, int num_bins = 5) {
    int n = x.size();
    if (n < 3) return NA_REAL;

    int lag = 2;
    int n_pairs = n - lag;

    double min_val = *std::min_element(x.begin(), x.end());
    double max_val = *std::max_element(x.begin(), x.end());

    if (max_val == min_val) return 0.0;

    double bin_width = (max_val - min_val) / num_bins;

    // Create 2D histogram for (x[t], x[t+lag])
    std::vector<std::vector<int>> joint_hist(num_bins, std::vector<int>(num_bins, 0));
    std::vector<int> marginal1(num_bins, 0);
    std::vector<int> marginal2(num_bins, 0);

    for (int i = 0; i < n_pairs; i++) {
        int bin1 = std::min(num_bins - 1, static_cast<int>((x[i] - min_val) / bin_width));
        int bin2 = std::min(num_bins - 1, static_cast<int>((x[i + lag] - min_val) / bin_width));

        joint_hist[bin1][bin2]++;
        marginal1[bin1]++;
        marginal2[bin2]++;
    }

    // Calculate mutual information
    double mi = 0.0;
    for (int i = 0; i < num_bins; i++) {
        for (int j = 0; j < num_bins; j++) {
            if (joint_hist[i][j] > 0) {
                double p_joint = static_cast<double>(joint_hist[i][j]) / n_pairs;
                double p_marg1 = static_cast<double>(marginal1[i]) / n_pairs;
                double p_marg2 = static_cast<double>(marginal2[j]) / n_pairs;

                mi += p_joint * std::log(p_joint / (p_marg1 * p_marg2));
            }
        }
    }

    return mi / std::log(2.0); // Convert to bits
}

// [[Rcpp::export]]
double cpp_automutual_info_first_min(NumericVector x, int max_lag = 40, int num_bins = 10) {
    int n = x.size();
    if (n < max_lag + 1) return NA_REAL;

    double min_val = *std::min_element(x.begin(), x.end());
    double max_val = *std::max_element(x.begin(), x.end());

    if (max_val == min_val) return NA_REAL;

    double bin_width = (max_val - min_val) / num_bins;

    std::vector<double> ami_values;

    for (int lag = 1; lag <= max_lag && lag < n; lag++) {
        int n_pairs = n - lag;

        // Create histograms
        std::vector<std::vector<int>> joint_hist(num_bins, std::vector<int>(num_bins, 0));
        std::vector<int> marginal1(num_bins, 0);
        std::vector<int> marginal2(num_bins, 0);

        for (int i = 0; i < n_pairs; i++) {
            int bin1 = std::min(num_bins - 1, static_cast<int>((x[i] - min_val) / bin_width));
            int bin2 = std::min(num_bins - 1, static_cast<int>((x[i + lag] - min_val) / bin_width));

            joint_hist[bin1][bin2]++;
            marginal1[bin1]++;
            marginal2[bin2]++;
        }

        // Calculate MI for this lag
        double mi = 0.0;
        for (int i = 0; i < num_bins; i++) {
            for (int j = 0; j < num_bins; j++) {
                if (joint_hist[i][j] > 0) {
                    double p_joint = static_cast<double>(joint_hist[i][j]) / n_pairs;
                    double p_marg1 = static_cast<double>(marginal1[i]) / n_pairs;
                    double p_marg2 = static_cast<double>(marginal2[j]) / n_pairs;

                    mi += p_joint * std::log(p_joint / (p_marg1 * p_marg2));
                }
            }
        }

        ami_values.push_back(mi);

        // Check for first minimum
        if (ami_values.size() >= 2) {
            int idx = ami_values.size() - 1;
            if (ami_values[idx] > ami_values[idx - 1] && idx >= 2) {
                return static_cast<double>(lag - 1);
            }
        }
    }

    return static_cast<double>(max_lag);
}

// [[Rcpp::export]]
double cpp_embedding_dist_exp_fit(NumericVector x) {
    int n = x.size();
    if (n < 10) return NA_REAL;

    // Simplified version: compute pairwise distances in 2D embedding
    // Using tau = 1 (lag-1 embedding)
    int m = 2; // embedding dimension
    int tau = 1;
    int n_vectors = n - (m - 1) * tau;

    if (n_vectors < 5) return NA_REAL;

    std::vector<double> distances;
    distances.reserve((n_vectors * (n_vectors - 1)) / 2);

    for (int i = 0; i < n_vectors; i++) {
        for (int j = i + 1; j < n_vectors; j++) {
            double dist = 0.0;
            for (int k = 0; k < m; k++) {
                double diff = x[i + k * tau] - x[j + k * tau];
                dist += diff * diff;
            }
            distances.push_back(std::sqrt(dist));
        }
    }

    std::sort(distances.begin(), distances.end());

    // Fit exponential distribution and return goodness of fit
    // Simplified: return mean distance / std distance as proxy
    double mean_dist = 0.0;
    for (double d : distances) mean_dist += d;
    mean_dist /= distances.size();

    double var_dist = 0.0;
    for (double d : distances) {
        var_dist += (d - mean_dist) * (d - mean_dist);
    }
    double std_dist = std::sqrt(var_dist / (distances.size() - 1));

    if (std_dist == 0.0) return NA_REAL;
    return mean_dist / std_dist; // Coefficient of variation (inverse)
}

// [[Rcpp::export]]
double cpp_dfa(NumericVector x) {
    int n = x.size();
    if (n < 20) return NA_REAL;

    // Integrate the signal (cumulative sum after removing mean)
    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= n;

    std::vector<double> y(n);
    y[0] = x[0] - mean;
    for (int i = 1; i < n; i++) {
        y[i] = y[i - 1] + (x[i] - mean);
    }

    // Divide into non-overlapping segments
    std::vector<int> box_sizes = {4, 8, 16, std::min(32, n / 4)};
    std::vector<double> log_box_sizes;
    std::vector<double> log_fluctuations;

    for (int box_size : box_sizes) {
        if (box_size >= n / 2) continue;

        int n_boxes = n / box_size;
        double sum_sq_residuals = 0.0;

        for (int box = 0; box < n_boxes; box++) {
            int start = box * box_size;
            int end = start + box_size;

            // Fit linear trend
            double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0;
            for (int i = start; i < end; i++) {
                double xi = i - start;
                sum_x += xi;
                sum_y += y[i];
                sum_xy += xi * y[i];
                sum_xx += xi * xi;
            }

            double slope = (box_size * sum_xy - sum_x * sum_y) / (box_size * sum_xx - sum_x * sum_x);
            double intercept = (sum_y - slope * sum_x) / box_size;

            // Calculate residuals
            for (int i = start; i < end; i++) {
                double xi = i - start;
                double trend = slope * xi + intercept;
                sum_sq_residuals += (y[i] - trend) * (y[i] - trend);
            }
        }

        double fluctuation = std::sqrt(sum_sq_residuals / (n_boxes * box_size));
        log_box_sizes.push_back(std::log(static_cast<double>(box_size)));
        log_fluctuations.push_back(std::log(fluctuation));
    }

    if (log_box_sizes.size() < 2) return NA_REAL;

    // Calculate slope (scaling exponent) using linear regression
    int m = log_box_sizes.size();
    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0;
    for (int i = 0; i < m; i++) {
        sum_x += log_box_sizes[i];
        sum_y += log_fluctuations[i];
        sum_xy += log_box_sizes[i] * log_fluctuations[i];
        sum_xx += log_box_sizes[i] * log_box_sizes[i];
    }

    double slope = (m * sum_xy - sum_x * sum_y) / (m * sum_xx - sum_x * sum_x);
    return slope;
}

// [[Rcpp::export]]
double cpp_rs_range(NumericVector x) {
    // Rescaled range analysis (R/S analysis)
    // Alternative to DFA for estimating Hurst exponent
    int n = x.size();
    if (n < 20) return NA_REAL;

    // Remove mean
    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= n;

    std::vector<double> y(n);
    for (int i = 0; i < n; i++) {
        y[i] = x[i] - mean;
    }

    // Calculate cumulative sum
    std::vector<double> cumsum(n);
    cumsum[0] = y[0];
    for (int i = 1; i < n; i++) {
        cumsum[i] = cumsum[i - 1] + y[i];
    }

    // Use different box sizes
    std::vector<int> box_sizes = {4, 8, 16, std::min(32, n / 4)};
    std::vector<double> log_box_sizes;
    std::vector<double> log_rs_values;

    for (int box_size : box_sizes) {
        if (box_size >= n / 2) continue;

        int n_boxes = n / box_size;
        double sum_rs = 0.0;

        for (int box = 0; box < n_boxes; box++) {
            int start = box * box_size;
            int end = start + box_size;

            // Calculate range of cumulative deviations within box
            double min_val = cumsum[start];
            double max_val = cumsum[start];

            for (int i = start + 1; i < end; i++) {
                if (cumsum[i] < min_val) min_val = cumsum[i];
                if (cumsum[i] > max_val) max_val = cumsum[i];
            }

            double range = max_val - min_val;

            // Calculate standard deviation within box
            double box_mean = 0.0;
            for (int i = start; i < end; i++) {
                box_mean += y[i];
            }
            box_mean /= box_size;

            double box_std = 0.0;
            for (int i = start; i < end; i++) {
                double diff = y[i] - box_mean;
                box_std += diff * diff;
            }
            box_std = std::sqrt(box_std / box_size);

            // R/S ratio (rescaled range)
            if (box_std > 0) {
                sum_rs += range / box_std;
            }
        }

        double avg_rs = sum_rs / n_boxes;
        log_box_sizes.push_back(std::log(static_cast<double>(box_size)));
        log_rs_values.push_back(std::log(avg_rs));
    }

    if (log_box_sizes.size() < 2) return NA_REAL;

    // Calculate Hurst exponent via linear regression
    // R/S ~ box_size^H, so log(R/S) ~ H * log(box_size)
    int m = log_box_sizes.size();
    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0;
    for (int i = 0; i < m; i++) {
        sum_x += log_box_sizes[i];
        sum_y += log_rs_values[i];
        sum_xy += log_box_sizes[i] * log_rs_values[i];
        sum_xx += log_box_sizes[i] * log_box_sizes[i];
    }

    double slope = (m * sum_xy - sum_x * sum_y) / (m * sum_xx - sum_x * sum_x);
    return slope; // Hurst exponent estimate
}

// [[Rcpp::export]]
double cpp_periodicity_wang(NumericVector x) {
    int n = x.size();
    if (n < 10) return NA_REAL;

    // Simplified Wang periodicity: based on ACF peaks
    // Calculate mean
    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= n;

    // Center the data
    std::vector<double> centered(n);
    for (int i = 0; i < n; i++) {
        centered[i] = x[i] - mean;
    }

    // Calculate variance
    double c0 = 0.0;
    for (int i = 0; i < n; i++) {
        c0 += centered[i] * centered[i];
    }

    if (c0 == 0.0) return 0.0;

    // Find ACF peaks
    int max_lag = std::min(n / 3, 100);
    std::vector<double> acf(max_lag);

    for (int lag = 1; lag < max_lag; lag++) {
        double c_lag = 0.0;
        for (int i = 0; i < n - lag; i++) {
            c_lag += centered[i] * centered[i + lag];
        }
        acf[lag] = c_lag / c0;
    }

    // Find highest ACF value (excluding lag 0)
    double max_acf = 0.0;
    for (int lag = 1; lag < max_lag; lag++) {
        if (acf[lag] > max_acf) {
            max_acf = acf[lag];
        }
    }

    return max_acf; // Simplified periodicity measure
}
