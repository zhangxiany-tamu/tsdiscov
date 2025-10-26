#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector cpp_acf(NumericVector x, int max_lag = 20, bool normalize = true) {
    int n = x.size();
    if (max_lag >= n) max_lag = n - 1;

    // Calculate mean
    double mean = 0.0;
    for (int i = 0; i < n; i++) {
        mean += x[i];
    }
    mean /= n;

    // Center the data
    NumericVector centered(n);
    for (int i = 0; i < n; i++) {
        centered[i] = x[i] - mean;
    }

    // Calculate autocorrelations
    NumericVector acf(max_lag + 1);
    double c0 = 0.0;

    // Variance (lag 0)
    for (int i = 0; i < n; i++) {
        c0 += centered[i] * centered[i];
    }
    acf[0] = 1.0;

    // Other lags
    for (int lag = 1; lag <= max_lag; lag++) {
        double c_lag = 0.0;
        for (int i = 0; i < n - lag; i++) {
            c_lag += centered[i] * centered[i + lag];
        }
        acf[lag] = normalize ? c_lag / c0 : c_lag;
    }

    return acf;
}

// [[Rcpp::export]]
double cpp_acf_first_min(NumericVector x, int max_lag = 20) {
    NumericVector acf = cpp_acf(x, max_lag, true);

    for (int i = 1; i < acf.size(); i++) {
        if (acf[i] < 0) {
            return static_cast<double>(i);
        }
    }
    return static_cast<double>(max_lag);
}

// [[Rcpp::export]]
double cpp_acf_timescale(NumericVector x, int max_lag = 20) {
    NumericVector acf = cpp_acf(x, max_lag, true);

    // Find where ACF first crosses exp(-1) ≈ 0.368
    double threshold = std::exp(-1.0);

    for (int i = 1; i < acf.size(); i++) {
        if (acf[i] < threshold) {
            // Linear interpolation
            double frac = (threshold - acf[i-1]) / (acf[i] - acf[i-1]);
            return (i - 1) + frac;
        }
    }
    return static_cast<double>(max_lag);
}

// [[Rcpp::export]]
NumericVector cpp_acf_features(NumericVector x, int max_lag = 20) {
    NumericVector acf = cpp_acf(x, max_lag, true);

    // Calculate various ACF-based features
    NumericVector features(6);

    // 1. First lag where ACF < 0
    features[0] = cpp_acf_first_min(x, max_lag);

    // 2. ACF timescale (crosses e^-1)
    features[1] = cpp_acf_timescale(x, max_lag);

    // 3. Sum of first 10 ACF values
    double sum_acf = 0.0;
    int limit = std::min(10, static_cast<int>(acf.size()) - 1);
    for (int i = 1; i <= limit; i++) {
        sum_acf += acf[i];
    }
    features[2] = sum_acf;

    // 4. ACF at lag 1
    features[3] = acf.size() > 1 ? acf[1] : NA_REAL;

    // 5. ACF at lag 10
    features[4] = acf.size() > 10 ? acf[10] : NA_REAL;

    // 6. Sum of squares of first 10 ACF values (x_acf10)
    double sum_sq_acf = 0.0;
    for (int i = 1; i <= limit; i++) {
        sum_sq_acf += acf[i] * acf[i];
    }
    features[5] = sum_sq_acf;

    return features;
}
