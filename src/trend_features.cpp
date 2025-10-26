#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector cpp_linear_trend(NumericVector x) {
    // Fit linear trend: y = slope * t + intercept
    // Returns: slope, intercept, R-squared, stderr
    int n = x.size();
    if (n < 3) return NumericVector::create(NA_REAL, NA_REAL, NA_REAL, NA_REAL);

    // Calculate sums for linear regression
    double sum_t = 0.0, sum_x = 0.0, sum_tx = 0.0, sum_tt = 0.0;
    for (int i = 0; i < n; i++) {
        double t = static_cast<double>(i);
        sum_t += t;
        sum_x += x[i];
        sum_tx += t * x[i];
        sum_tt += t * t;
    }

    // Calculate slope and intercept
    double n_dbl = static_cast<double>(n);
    double slope = (n_dbl * sum_tx - sum_t * sum_x) / (n_dbl * sum_tt - sum_t * sum_t);
    double intercept = (sum_x - slope * sum_t) / n_dbl;

    // Calculate R-squared
    double mean_x = sum_x / n_dbl;
    double ss_tot = 0.0, ss_res = 0.0;
    for (int i = 0; i < n; i++) {
        double t = static_cast<double>(i);
        double y_pred = slope * t + intercept;
        ss_tot += (x[i] - mean_x) * (x[i] - mean_x);
        ss_res += (x[i] - y_pred) * (x[i] - y_pred);
    }

    double r_squared = (ss_tot > 0) ? 1.0 - (ss_res / ss_tot) : NA_REAL;

    // Calculate standard error of slope
    double stderr_slope = NA_REAL;
    if (n > 2 && ss_tot > 0) {
        double mse = ss_res / (n - 2);
        double s_tt = sum_tt - (sum_t * sum_t) / n_dbl;
        stderr_slope = std::sqrt(mse / s_tt);
    }

    return NumericVector::create(
        Named("slope") = slope,
        Named("intercept") = intercept,
        Named("r_squared") = r_squared,
        Named("stderr") = stderr_slope
    );
}

// [[Rcpp::export]]
double cpp_mean_abs_change(NumericVector x) {
    // Mean of absolute differences between consecutive values
    int n = x.size();
    if (n < 2) return NA_REAL;

    double sum_abs_diff = 0.0;
    for (int i = 1; i < n; i++) {
        sum_abs_diff += std::abs(x[i] - x[i - 1]);
    }

    return sum_abs_diff / (n - 1);
}

// [[Rcpp::export]]
double cpp_mean_change(NumericVector x) {
    // Mean of differences between consecutive values
    int n = x.size();
    if (n < 2) return NA_REAL;

    double sum_diff = 0.0;
    for (int i = 1; i < n; i++) {
        sum_diff += x[i] - x[i - 1];
    }

    return sum_diff / (n - 1);
}

// [[Rcpp::export]]
double cpp_mean_second_derivative(NumericVector x) {
    // Mean of second derivative (acceleration)
    int n = x.size();
    if (n < 3) return NA_REAL;

    double sum_second_deriv = 0.0;
    for (int i = 1; i < n - 1; i++) {
        // Second derivative: (x[i+1] - x[i]) - (x[i] - x[i-1])
        sum_second_deriv += (x[i + 1] - 2.0 * x[i] + x[i - 1]) / 2.0;
    }

    return sum_second_deriv / (n - 2);
}

// Time-weighted linear trend - recent data weighted more heavily
// [[Rcpp::export]]
NumericVector cpp_time_weighted_trend(NumericVector x, double decay = 0.95) {
    int n = x.size();

    // Default result
    NumericVector default_result = NumericVector::create(
        Named("time_weighted_slope") = NA_REAL,
        Named("time_weighted_intercept") = NA_REAL
    );

    if (n < 10) return default_result;

    // Create exponential weights (recent data gets higher weight)
    NumericVector weights(n);
    double weight_sum = 0.0;
    for (int i = 0; i < n; i++) {
        // Weight increases exponentially toward the end
        weights[i] = std::pow(decay, n - 1 - i);
        weight_sum += weights[i];
    }

    // Normalize weights
    for (int i = 0; i < n; i++) {
        weights[i] /= weight_sum;
    }

    // Time index
    NumericVector time(n);
    for (int i = 0; i < n; i++) {
        time[i] = static_cast<double>(i);
    }

    // Weighted means
    double mean_x = 0.0, mean_t = 0.0;
    for (int i = 0; i < n; i++) {
        mean_x += weights[i] * x[i];
        mean_t += weights[i] * time[i];
    }

    // Weighted covariance and variance
    double cov = 0.0, var_t = 0.0;
    for (int i = 0; i < n; i++) {
        cov += weights[i] * (time[i] - mean_t) * (x[i] - mean_x);
        var_t += weights[i] * (time[i] - mean_t) * (time[i] - mean_t);
    }

    if (var_t < 1e-10) return default_result;

    double slope = cov / var_t;
    double intercept = mean_x - slope * mean_t;

    return NumericVector::create(
        Named("time_weighted_slope") = slope,
        Named("time_weighted_intercept") = intercept
    );
}

// Theil-Sen robust trend estimator - median of all pairwise slopes
// [[Rcpp::export]]
NumericVector cpp_robust_trend(NumericVector x, int max_pairs = 5000) {
    int n = x.size();

    // Default result
    NumericVector default_result = NumericVector::create(
        Named("robust_slope") = NA_REAL,
        Named("robust_intercept") = NA_REAL,
        Named("robust_slope_std") = NA_REAL
    );

    if (n < 10) return default_result;

    // Create time indices
    std::vector<int> time_idx(n);
    for (int i = 0; i < n; i++) {
        time_idx[i] = i;
    }

    // Calculate slopes for all pairs (or sample if too many)
    std::vector<double> slopes;

    int total_pairs = n * (n - 1) / 2;
    bool sample = (total_pairs > max_pairs);

    if (sample) {
        // Sample pairs randomly
        slopes.reserve(max_pairs);
        std::srand(123);  // Fixed seed for reproducibility

        for (int k = 0; k < max_pairs; k++) {
            int i = std::rand() % n;
            int j = std::rand() % n;
            if (i == j) continue;

            if (time_idx[j] != time_idx[i]) {  // Avoid division by zero
                double slope = (x[j] - x[i]) / (time_idx[j] - time_idx[i]);
                slopes.push_back(slope);
            }
        }
    } else {
        // Calculate all pairwise slopes
        slopes.reserve(total_pairs);

        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                double slope = (x[j] - x[i]) / (time_idx[j] - time_idx[i]);
                slopes.push_back(slope);
            }
        }
    }

    if (slopes.empty()) return default_result;

    // Compute median slope (Theil-Sen estimator)
    std::sort(slopes.begin(), slopes.end());
    double median_slope;
    int n_slopes = slopes.size();

    if (n_slopes % 2 == 0) {
        median_slope = (slopes[n_slopes / 2 - 1] + slopes[n_slopes / 2]) / 2.0;
    } else {
        median_slope = slopes[n_slopes / 2];
    }

    // Compute intercept: median(y - slope * x)
    std::vector<double> intercepts(n);
    for (int i = 0; i < n; i++) {
        intercepts[i] = x[i] - median_slope * time_idx[i];
    }

    std::sort(intercepts.begin(), intercepts.end());
    double median_intercept;

    if (n % 2 == 0) {
        median_intercept = (intercepts[n / 2 - 1] + intercepts[n / 2]) / 2.0;
    } else {
        median_intercept = intercepts[n / 2];
    }

    // Compute standard deviation of slopes (measure of uncertainty)
    double mean_slope = 0.0;
    for (double s : slopes) {
        mean_slope += s;
    }
    mean_slope /= n_slopes;

    double var_slope = 0.0;
    for (double s : slopes) {
        var_slope += (s - mean_slope) * (s - mean_slope);
    }
    double std_slope = std::sqrt(var_slope / n_slopes);

    return NumericVector::create(
        Named("robust_slope") = median_slope,
        Named("robust_intercept") = median_intercept,
        Named("robust_slope_std") = std_slope
    );
}
