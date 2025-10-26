#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
double cpp_sample_entropy(NumericVector x, int m = 2, double r = 0.2) {
    int n = x.size();
    if (n < m + 1) return NA_REAL;

    // Calculate standard deviation for threshold
    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= n;

    double sd = 0.0;
    for (int i = 0; i < n; i++) {
        sd += (x[i] - mean) * (x[i] - mean);
    }
    sd = std::sqrt(sd / (n - 1));

    double threshold = r * sd;

    // Count template matches
    auto count_matches = [&](int template_length) -> double {
        int count = 0;
        int N = n - template_length;

        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                bool match = true;
                for (int k = 0; k < template_length; k++) {
                    if (std::abs(x[i + k] - x[j + k]) > threshold) {
                        match = false;
                        break;
                    }
                }
                if (match) count++;
            }
        }
        return static_cast<double>(count);
    };

    double A = count_matches(m);
    double B = count_matches(m + 1);

    if (A == 0 || B == 0) return NA_REAL;

    return -std::log(B / A);
}

// [[Rcpp::export]]
double cpp_approximate_entropy(NumericVector x, int m = 2, double r = 0.2) {
    int n = x.size();
    if (n < m + 1) return NA_REAL;

    // Calculate standard deviation for threshold
    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= n;

    double sd = 0.0;
    for (int i = 0; i < n; i++) {
        sd += (x[i] - mean) * (x[i] - mean);
    }
    sd = std::sqrt(sd / (n - 1));

    double threshold = r * sd;

    auto calculate_phi = [&](int template_length) -> double {
        std::vector<double> C(n - template_length + 1, 0.0);

        for (int i = 0; i <= n - template_length; i++) {
            int count = 0;
            for (int j = 0; j <= n - template_length; j++) {
                bool match = true;
                for (int k = 0; k < template_length; k++) {
                    if (std::abs(x[i + k] - x[j + k]) > threshold) {
                        match = false;
                        break;
                    }
                }
                if (match) count++;
            }
            C[i] = static_cast<double>(count) / (n - template_length + 1);
        }

        double phi = 0.0;
        for (size_t i = 0; i < C.size(); i++) {
            if (C[i] > 0) phi += std::log(C[i]);
        }
        return phi / C.size();
    };

    double phi_m = calculate_phi(m);
    double phi_m1 = calculate_phi(m + 1);

    return phi_m - phi_m1;
}

// [[Rcpp::export]]
double cpp_permutation_entropy(NumericVector x, int m = 3, int tau = 1) {
    int n = x.size();
    int num_patterns = n - (m - 1) * tau;

    if (num_patterns < 1) return NA_REAL;

    std::map<std::vector<int>, int> pattern_counts;

    // Extract patterns
    for (int i = 0; i < num_patterns; i++) {
        std::vector<double> values(m);
        std::vector<int> indices(m);

        for (int j = 0; j < m; j++) {
            values[j] = x[i + j * tau];
            indices[j] = j;
        }

        // Sort indices by values
        std::sort(indices.begin(), indices.end(),
                  [&values](int a, int b) { return values[a] < values[b]; });

        pattern_counts[indices]++;
    }

    // Calculate entropy
    double entropy = 0.0;
    for (const auto& pair : pattern_counts) {
        double p = static_cast<double>(pair.second) / num_patterns;
        entropy -= p * std::log(p);
    }

    return entropy / std::log(2.0);  // Normalize to bits
}

// [[Rcpp::export]]
double cpp_shannon_entropy(NumericVector x, int num_bins = 10) {
    int n = x.size();
    if (n == 0) return NA_REAL;

    double min_val = *std::min_element(x.begin(), x.end());
    double max_val = *std::max_element(x.begin(), x.end());

    if (max_val == min_val) return 0.0;

    std::vector<int> bins(num_bins, 0);
    double bin_width = (max_val - min_val) / num_bins;

    // Fill histogram
    for (int i = 0; i < n; i++) {
        int bin = static_cast<int>((x[i] - min_val) / bin_width);
        if (bin >= num_bins) bin = num_bins - 1;
        bins[bin]++;
    }

    // Calculate entropy
    double entropy = 0.0;
    for (int i = 0; i < num_bins; i++) {
        if (bins[i] > 0) {
            double p = static_cast<double>(bins[i]) / n;
            entropy -= p * std::log(p);
        }
    }

    return entropy / std::log(2.0);  // Normalize to bits
}

// Permutation Entropy (Bandt-Pompe) - measures complexity via ordinal patterns
// [[Rcpp::export]]
NumericVector cpp_multiscale_entropy(NumericVector x, IntegerVector scales = IntegerVector::create(2, 5), int m = 2, double r = 0.2) {
    int n = x.size();
    int n_scales = scales.size();

    NumericVector result(n_scales);
    CharacterVector names(n_scales);

    for (int s_idx = 0; s_idx < n_scales; s_idx++) {
        int scale = scales[s_idx];
        names[s_idx] = "multiscale_entropy_" + std::to_string(scale);

        // Coarse-grain at this scale
        int n_coarse = n / scale;
        if (n_coarse < m + 10) {
            result[s_idx] = NA_REAL;
            continue;
        }

        NumericVector coarse(n_coarse);
        for (int i = 0; i < n_coarse; i++) {
            double sum = 0.0;
            for (int j = 0; j < scale; j++) {
                sum += x[i * scale + j];
            }
            coarse[i] = sum / scale;
        }

        // Compute sample entropy on coarse-grained series
        result[s_idx] = cpp_sample_entropy(coarse, m, r);
    }

    result.names() = names;
    return result;
}
