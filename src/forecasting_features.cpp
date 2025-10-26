#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
double cpp_forecast_error_mean3(NumericVector x) {
    int n = x.size();
    if (n < 6) return NA_REAL;

    // Simple 3-point rolling mean forecast
    // Use first half to forecast second half
    int train_size = n / 2;
    double sum_sq_error = 0.0;
    int count = 0;

    for (int i = train_size; i < n; i++) {
        if (i >= 3) {
            // Forecast as mean of previous 3 points
            double forecast = (x[i-1] + x[i-2] + x[i-3]) / 3.0;
            double error = x[i] - forecast;
            sum_sq_error += error * error;
            count++;
        }
    }

    if (count == 0) return NA_REAL;
    return std::sqrt(sum_sq_error / count); // RMSE
}

// Helper: compute ACF timescale
double compute_acf_timescale(NumericVector x, int max_lag) {
    int n = x.size();

    // Calculate mean
    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= n;

    // Center the data
    NumericVector centered(n);
    for (int i = 0; i < n; i++) {
        centered[i] = x[i] - mean;
    }

    // Calculate variance (lag 0)
    double c0 = 0.0;
    for (int i = 0; i < n; i++) {
        c0 += centered[i] * centered[i];
    }

    if (c0 == 0.0) return NA_REAL;

    // Find where ACF crosses exp(-1) ≈ 0.368
    double threshold = std::exp(-1.0);

    for (int lag = 1; lag <= max_lag && lag < n; lag++) {
        double c_lag = 0.0;
        for (int i = 0; i < n - lag; i++) {
            c_lag += centered[i] * centered[i + lag];
        }
        double acf = c_lag / c0;

        if (acf < threshold) {
            // Linear interpolation
            double c_prev = 0.0;
            for (int i = 0; i < n - (lag - 1); i++) {
                c_prev += centered[i] * centered[i + lag - 1];
            }
            double acf_prev = c_prev / c0;

            double frac = (threshold - acf_prev) / (acf - acf_prev);
            return (lag - 1) + frac;
        }
    }

    return static_cast<double>(max_lag);
}

// [[Rcpp::export]]
double cpp_timescale_ratio_after_whitening(NumericVector x) {
    int n = x.size();
    if (n < 4) return NA_REAL;

    // Calculate original timescale
    double tau_original = compute_acf_timescale(x, std::min(20, n - 1));

    // Difference the series (simple whitening)
    NumericVector diff(n - 1);
    for (int i = 0; i < n - 1; i++) {
        diff[i] = x[i + 1] - x[i];
    }

    // Calculate differenced timescale
    double tau_diff = compute_acf_timescale(diff, std::min(20, static_cast<int>(diff.size()) - 1));

    if (tau_original == 0.0) return NA_REAL;
    return tau_diff / tau_original;
}

// [[Rcpp::export]]
double cpp_time_reversibility(NumericVector x) {
    int n = x.size();
    if (n < 3) return NA_REAL;

    // Time reversibility statistic: E[(x[t+1] - x[t])^3]
    // Measures asymmetry under time reversal
    double sum_cube = 0.0;
    int count = 0;

    for (int i = 0; i < n - 1; i++) {
        double diff = x[i + 1] - x[i];
        sum_cube += diff * diff * diff;
        count++;
    }

    if (count == 0) return NA_REAL;
    return sum_cube / count;
}

// [[Rcpp::export]]
double cpp_high_fluctuation_prop(NumericVector x) {
    // MD_hrv_classic_pnn40: proportion of incremental changes > 40% of std
    int n = x.size();
    if (n < 2) return NA_REAL;

    // Calculate standard deviation of differences
    NumericVector diff(n - 1);
    for (int i = 0; i < n - 1; i++) {
        diff[i] = x[i + 1] - x[i];
    }

    double mean_diff = 0.0;
    for (int i = 0; i < n - 1; i++) {
        mean_diff += diff[i];
    }
    mean_diff /= (n - 1);

    double var_diff = 0.0;
    for (int i = 0; i < n - 1; i++) {
        double d = diff[i] - mean_diff;
        var_diff += d * d;
    }
    double std_diff = std::sqrt(var_diff / (n - 2));

    if (std_diff == 0.0) return 0.0;

    // Count proportion of changes > 0.4 * std
    double threshold = 0.4 * std_diff;
    int count = 0;
    for (int i = 0; i < n - 1; i++) {
        if (std::abs(diff[i]) > threshold) {
            count++;
        }
    }

    return static_cast<double>(count) / (n - 1);
}
