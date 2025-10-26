#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// Helper function: Compute ACF
std::vector<double> compute_acf_vector(const NumericVector& x, int max_lag) {
    int n = x.size();
    if (max_lag >= n) max_lag = n - 1;

    // Calculate mean
    double mean = 0.0;
    for (int i = 0; i < n; i++) {
        mean += x[i];
    }
    mean /= n;

    // Center the data
    std::vector<double> centered(n);
    for (int i = 0; i < n; i++) {
        centered[i] = x[i] - mean;
    }

    // Calculate variance (lag 0)
    double c0 = 0.0;
    for (int i = 0; i < n; i++) {
        c0 += centered[i] * centered[i];
    }

    if (c0 == 0.0) return std::vector<double>(max_lag + 1, 0.0);

    // Calculate ACF
    std::vector<double> acf(max_lag + 1);
    acf[0] = 1.0;

    for (int lag = 1; lag <= max_lag; lag++) {
        double c_lag = 0.0;
        for (int i = 0; i < n - lag; i++) {
            c_lag += centered[i] * centered[i + lag];
        }
        acf[lag] = c_lag / c0;
    }

    return acf;
}

// [[Rcpp::export]]
NumericVector cpp_pacf(NumericVector x, int max_lag = 20) {
    // Compute partial autocorrelation function using Durbin-Levinson algorithm
    int n = x.size();
    if (n < 3 || max_lag < 1) return NumericVector::create(NA_REAL);
    if (max_lag >= n) max_lag = n - 1;

    // Get ACF values
    std::vector<double> acf = compute_acf_vector(x, max_lag);

    // Durbin-Levinson algorithm for PACF
    std::vector<double> pacf(max_lag + 1);
    pacf[0] = 1.0;

    if (max_lag >= 1) {
        pacf[1] = acf[1];
    }

    std::vector<std::vector<double>> phi(max_lag + 1);
    for (int i = 0; i <= max_lag; i++) {
        phi[i].resize(i + 1);
    }

    phi[1][1] = acf[1];

    for (int k = 2; k <= max_lag; k++) {
        // Calculate partial autocorrelation at lag k
        double numerator = acf[k];
        double denominator = 1.0;

        for (int j = 1; j < k; j++) {
            numerator -= phi[k - 1][j] * acf[k - j];
            denominator -= phi[k - 1][j] * acf[j];
        }

        if (std::abs(denominator) < 1e-10) {
            pacf[k] = 0.0;
        } else {
            pacf[k] = numerator / denominator;
        }

        // Update phi coefficients
        phi[k][k] = pacf[k];
        for (int j = 1; j < k; j++) {
            phi[k][j] = phi[k - 1][j] - pacf[k] * phi[k - 1][k - j];
        }
    }

    // Return PACF values (excluding lag 0)
    NumericVector result(max_lag);
    for (int i = 1; i <= max_lag; i++) {
        result[i - 1] = pacf[i];
    }

    return result;
}

// [[Rcpp::export]]
NumericVector cpp_pacf_features(NumericVector x, int max_lag = 20) {
    // Extract key PACF-based features
    NumericVector pacf = cpp_pacf(x, max_lag);

    if (pacf.size() < 10) {
        return NumericVector::create(NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL);
    }

    // Feature 1: PACF at lag 1
    double pacf_lag1 = pacf[0];

    // Feature 2: PACF at lag 5
    double pacf_lag5 = (pacf.size() >= 5) ? pacf[4] : NA_REAL;

    // Feature 3: PACF at lag 10
    double pacf_lag10 = (pacf.size() >= 10) ? pacf[9] : NA_REAL;

    // Feature 4: First significant PACF lag (|PACF| > 2/sqrt(n))
    int n = x.size();
    double threshold = 2.0 / std::sqrt(static_cast<double>(n));
    double first_sig_lag = 0.0;

    for (int i = 0; i < pacf.size(); i++) {
        if (std::abs(pacf[i]) > threshold) {
            first_sig_lag = static_cast<double>(i + 1);
            break;
        }
    }

    if (first_sig_lag == 0.0) {
        first_sig_lag = static_cast<double>(max_lag);
    }

    // Feature 5: Sum of squared PACF values (first 10)
    double sum_sq_pacf = 0.0;
    int limit = std::min(10, static_cast<int>(pacf.size()));
    for (int i = 0; i < limit; i++) {
        sum_sq_pacf += pacf[i] * pacf[i];
    }

    // Feature 6: Sum of squared PACF values (first 5) - x_pacf5
    double x_pacf5 = 0.0;
    int limit5 = std::min(5, static_cast<int>(pacf.size()));
    for (int i = 0; i < limit5; i++) {
        x_pacf5 += pacf[i] * pacf[i];
    }

    return NumericVector::create(
        Named("pacf_lag1") = pacf_lag1,
        Named("pacf_lag5") = pacf_lag5,
        Named("pacf_lag10") = pacf_lag10,
        Named("first_sig_pacf") = first_sig_lag,
        Named("sum_sq_pacf") = sum_sq_pacf,
        Named("x_pacf5") = x_pacf5
    );
}
