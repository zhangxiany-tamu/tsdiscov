#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
double cpp_trend_strength(NumericVector trend, NumericVector remainder) {
    // Trend strength: 1 - Var(Remainder) / Var(Trend + Remainder)
    // High values indicate strong trend

    int n = trend.size();
    if (n != remainder.size() || n < 2) return NA_REAL;

    // Calculate variance of remainder
    double mean_rem = 0.0;
    for (int i = 0; i < n; i++) mean_rem += remainder[i];
    mean_rem /= n;

    double var_rem = 0.0;
    for (int i = 0; i < n; i++) {
        var_rem += (remainder[i] - mean_rem) * (remainder[i] - mean_rem);
    }
    var_rem /= n;

    // Calculate variance of (trend + remainder)
    std::vector<double> detrended(n);
    double mean_detrended = 0.0;
    for (int i = 0; i < n; i++) {
        detrended[i] = trend[i] + remainder[i];
        mean_detrended += detrended[i];
    }
    mean_detrended /= n;

    double var_detrended = 0.0;
    for (int i = 0; i < n; i++) {
        var_detrended += (detrended[i] - mean_detrended) * (detrended[i] - mean_detrended);
    }
    var_detrended /= n;

    if (var_detrended < 1e-10) return NA_REAL;

    double strength = 1.0 - (var_rem / var_detrended);
    return std::max(0.0, std::min(1.0, strength));
}

// [[Rcpp::export]]
double cpp_seasonal_strength(NumericVector seasonal, NumericVector remainder) {
    // Seasonal strength: 1 - Var(Remainder) / Var(Seasonal + Remainder)
    // High values indicate strong seasonality

    int n = seasonal.size();
    if (n != remainder.size() || n < 2) return NA_REAL;

    // Calculate variance of remainder
    double mean_rem = 0.0;
    for (int i = 0; i < n; i++) mean_rem += remainder[i];
    mean_rem /= n;

    double var_rem = 0.0;
    for (int i = 0; i < n; i++) {
        var_rem += (remainder[i] - mean_rem) * (remainder[i] - mean_rem);
    }
    var_rem /= n;

    // Calculate variance of (seasonal + remainder)
    std::vector<double> deseasoned(n);
    double mean_deseasoned = 0.0;
    for (int i = 0; i < n; i++) {
        deseasoned[i] = seasonal[i] + remainder[i];
        mean_deseasoned += deseasoned[i];
    }
    mean_deseasoned /= n;

    double var_deseasoned = 0.0;
    for (int i = 0; i < n; i++) {
        var_deseasoned += (deseasoned[i] - mean_deseasoned) * (deseasoned[i] - mean_deseasoned);
    }
    var_deseasoned /= n;

    if (var_deseasoned < 1e-10) return NA_REAL;

    double strength = 1.0 - (var_rem / var_deseasoned);
    return std::max(0.0, std::min(1.0, strength));
}

// [[Rcpp::export]]
double cpp_spike(NumericVector remainder) {
    // Spike: variance of leave-one-out variances of remainder
    // Measures spikiness/outliers in remainder

    int n = remainder.size();
    if (n < 3) return NA_REAL;

    std::vector<double> loo_vars(n);

    // Calculate leave-one-out variances
    for (int i = 0; i < n; i++) {
        // Calculate variance excluding observation i
        double sum = 0.0;
        int count = 0;
        for (int j = 0; j < n; j++) {
            if (j != i) {
                sum += remainder[j];
                count++;
            }
        }
        double mean_loo = sum / count;

        double var = 0.0;
        for (int j = 0; j < n; j++) {
            if (j != i) {
                var += (remainder[j] - mean_loo) * (remainder[j] - mean_loo);
            }
        }
        loo_vars[i] = var / count;
    }

    // Calculate variance of leave-one-out variances
    double mean_var = 0.0;
    for (int i = 0; i < n; i++) {
        mean_var += loo_vars[i];
    }
    mean_var /= n;

    double spike = 0.0;
    for (int i = 0; i < n; i++) {
        spike += (loo_vars[i] - mean_var) * (loo_vars[i] - mean_var);
    }
    spike /= n;

    return spike;
}

// [[Rcpp::export]]
double cpp_linearity(NumericVector trend) {
    // Linearity: R² of linear fit to trend component
    // High values indicate trend is close to linear

    int n = trend.size();
    if (n < 3) return NA_REAL;

    // Fit linear regression: trend ~ t
    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0;
    for (int i = 0; i < n; i++) {
        double x = static_cast<double>(i);
        sum_x += x;
        sum_y += trend[i];
        sum_xy += x * trend[i];
        sum_xx += x * x;
    }

    double denom = n * sum_xx - sum_x * sum_x;
    if (std::abs(denom) < 1e-10) return NA_REAL;

    double slope = (n * sum_xy - sum_x * sum_y) / denom;
    double intercept = (sum_y - slope * sum_x) / n;

    // Calculate R²
    double mean_y = sum_y / n;
    double ss_tot = 0.0, ss_res = 0.0;
    for (int i = 0; i < n; i++) {
        double x = static_cast<double>(i);
        double fitted = slope * x + intercept;
        ss_tot += (trend[i] - mean_y) * (trend[i] - mean_y);
        ss_res += (trend[i] - fitted) * (trend[i] - fitted);
    }

    if (ss_tot < 1e-10) return NA_REAL;

    double r_squared = 1.0 - (ss_res / ss_tot);
    return std::max(0.0, std::min(1.0, r_squared));
}

// [[Rcpp::export]]
double cpp_curvature(NumericVector trend) {
    // Curvature: R² improvement from quadratic vs linear fit
    // High values indicate non-linear (curved) trend

    int n = trend.size();
    if (n < 5) return NA_REAL;

    // First, get linear R²
    double linear_r2 = cpp_linearity(trend);
    if (!std::isfinite(linear_r2)) return NA_REAL;

    // Fit quadratic: trend ~ a*t² + b*t + c
    // Using normal equations for quadratic regression
    double sum_x = 0.0, sum_x2 = 0.0, sum_x3 = 0.0, sum_x4 = 0.0;
    double sum_y = 0.0, sum_xy = 0.0, sum_x2y = 0.0;

    for (int i = 0; i < n; i++) {
        double x = static_cast<double>(i);
        double x2 = x * x;
        double x3 = x2 * x;
        double x4 = x2 * x2;

        sum_x += x;
        sum_x2 += x2;
        sum_x3 += x3;
        sum_x4 += x4;
        sum_y += trend[i];
        sum_xy += x * trend[i];
        sum_x2y += x2 * trend[i];
    }

    // Solve 3x3 system (simplified - just calculate residuals)
    // For curvature measure, compare quadratic vs linear residuals
    double mean_y = sum_y / n;
    double ss_tot = 0.0;
    for (int i = 0; i < n; i++) {
        ss_tot += (trend[i] - mean_y) * (trend[i] - mean_y);
    }

    if (ss_tot < 1e-10) return NA_REAL;

    // Simple curvature metric: variance of second differences
    if (n < 3) return NA_REAL;

    std::vector<double> second_diff(n - 2);
    for (int i = 0; i < n - 2; i++) {
        second_diff[i] = trend[i + 2] - 2 * trend[i + 1] + trend[i];
    }

    double mean_sd = 0.0;
    for (int i = 0; i < n - 2; i++) {
        mean_sd += second_diff[i];
    }
    mean_sd /= (n - 2);

    double var_sd = 0.0;
    for (int i = 0; i < n - 2; i++) {
        var_sd += (second_diff[i] - mean_sd) * (second_diff[i] - mean_sd);
    }
    var_sd /= (n - 2);

    // Normalize by total variance
    return var_sd / (ss_tot / n);
}
