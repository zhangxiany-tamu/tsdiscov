#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// Forward declaration of cpp_rs_range from advanced_features.cpp
double cpp_rs_range(NumericVector x);

// [[Rcpp::export]]
double cpp_hurst_exponent(NumericVector x) {
    // Calculate Hurst exponent using R/S range method
    // This is already implemented in advanced_features.cpp as cpp_rs_range
    // H ~ 0.5: white noise
    // H > 0.5: persistent (trending)
    // H < 0.5: anti-persistent (mean-reverting)

    return cpp_rs_range(x);
}

// [[Rcpp::export]]
double cpp_stability(NumericVector x, int window_size = 10) {
    // Stability: variance of means in non-overlapping windows
    // Low stability = high variance in local means (non-stationary)
    // High stability = low variance in local means (stationary)

    int n = x.size();
    if (n < window_size * 2) return NA_REAL;

    int n_windows = n / window_size;
    if (n_windows < 2) return NA_REAL;

    std::vector<double> window_means(n_windows);

    // Calculate mean for each window
    for (int i = 0; i < n_windows; i++) {
        double sum = 0.0;
        for (int j = 0; j < window_size; j++) {
            sum += x[i * window_size + j];
        }
        window_means[i] = sum / window_size;
    }

    // Calculate variance of window means
    double mean_of_means = 0.0;
    for (int i = 0; i < n_windows; i++) {
        mean_of_means += window_means[i];
    }
    mean_of_means /= n_windows;

    double variance = 0.0;
    for (int i = 0; i < n_windows; i++) {
        variance += (window_means[i] - mean_of_means) * (window_means[i] - mean_of_means);
    }
    variance /= n_windows;

    return variance;
}

// [[Rcpp::export]]
double cpp_lumpiness(NumericVector x, int window_size = 10) {
    // Lumpiness: variance of variances in non-overlapping windows
    // High lumpiness = variance changes over time (heteroscedastic)
    // Low lumpiness = stable variance (homoscedastic)

    int n = x.size();
    if (n < window_size * 2) return NA_REAL;

    int n_windows = n / window_size;
    if (n_windows < 2) return NA_REAL;

    std::vector<double> window_vars(n_windows);

    // Calculate variance for each window
    for (int i = 0; i < n_windows; i++) {
        // First calculate mean of window
        double window_mean = 0.0;
        for (int j = 0; j < window_size; j++) {
            window_mean += x[i * window_size + j];
        }
        window_mean /= window_size;

        // Then calculate variance
        double var = 0.0;
        for (int j = 0; j < window_size; j++) {
            double diff = x[i * window_size + j] - window_mean;
            var += diff * diff;
        }
        window_vars[i] = var / window_size;
    }

    // Calculate variance of window variances
    double mean_of_vars = 0.0;
    for (int i = 0; i < n_windows; i++) {
        mean_of_vars += window_vars[i];
    }
    mean_of_vars /= n_windows;

    double variance = 0.0;
    for (int i = 0; i < n_windows; i++) {
        variance += (window_vars[i] - mean_of_vars) * (window_vars[i] - mean_of_vars);
    }
    variance /= n_windows;

    return variance;
}
