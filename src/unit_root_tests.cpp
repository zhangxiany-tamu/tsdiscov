#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// Helper function for OLS regression
// Returns coefficients and standard errors
List ols_regression(NumericMatrix X, NumericVector y) {
    int n = X.nrow();
    int k = X.ncol();

    // Check for sufficient data
    if (n <= k) {
        return List::create(
            Named("coefficients") = NumericVector(k, NA_REAL),
            Named("std_errors") = NumericVector(k, NA_REAL),
            Named("residuals") = NumericVector(n, NA_REAL)
        );
    }

    // Compute X'X
    NumericMatrix XtX(k, k);
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            double sum = 0.0;
            for (int row = 0; row < n; row++) {
                sum += X(row, i) * X(row, j);
            }
            XtX(i, j) = sum;
        }
    }

    // Compute X'y
    NumericVector Xty(k);
    for (int i = 0; i < k; i++) {
        double sum = 0.0;
        for (int row = 0; row < n; row++) {
            sum += X(row, i) * y[row];
        }
        Xty[i] = sum;
    }

    // Solve (X'X)^{-1} X'y using Cholesky decomposition
    // For numerical stability, use R's built-in solve
    Environment base("package:base");
    Function solve_R = base["solve"];

    NumericVector coefficients;
    try {
        coefficients = as<NumericVector>(solve_R(XtX, Xty));
    } catch(...) {
        // If matrix is singular, return NAs
        return List::create(
            Named("coefficients") = NumericVector(k, NA_REAL),
            Named("std_errors") = NumericVector(k, NA_REAL),
            Named("residuals") = NumericVector(n, NA_REAL)
        );
    }

    // Compute residuals
    NumericVector residuals(n);
    for (int i = 0; i < n; i++) {
        double fitted = 0.0;
        for (int j = 0; j < k; j++) {
            fitted += X(i, j) * coefficients[j];
        }
        residuals[i] = y[i] - fitted;
    }

    // Compute residual sum of squares
    double rss = 0.0;
    for (int i = 0; i < n; i++) {
        rss += residuals[i] * residuals[i];
    }

    // Estimate error variance
    double sigma2 = rss / (n - k);

    // Compute standard errors using (X'X)^{-1} * sigma2
    NumericMatrix XtX_inv;
    try {
        XtX_inv = as<NumericMatrix>(solve_R(XtX));
    } catch(...) {
        return List::create(
            Named("coefficients") = coefficients,
            Named("std_errors") = NumericVector(k, NA_REAL),
            Named("residuals") = residuals
        );
    }

    NumericVector std_errors(k);
    for (int i = 0; i < k; i++) {
        std_errors[i] = std::sqrt(sigma2 * XtX_inv(i, i));
    }

    return List::create(
        Named("coefficients") = coefficients,
        Named("std_errors") = std_errors,
        Named("residuals") = residuals
    );
}

// MacKinnon (1994) p-value approximation for ADF test
// For case with constant and trend (tau_ct)
double adf_pvalue_mackinnon(double test_stat, int n) {
    // MacKinnon's coefficients for tau_ct (constant + trend)
    // For T -> infinity: critical values at 1%, 5%, 10%: -3.96, -3.41, -3.13

    // Small sample correction factors
    double beta_inf, beta_1, beta_2;

    // Asymptotic coefficients for tau_ct
    beta_inf = -3.4326;
    beta_1 = -5.999;
    beta_2 = -29.25;

    // Small sample adjustment
    double cv = beta_inf + beta_1 / n + beta_2 / (n * n);

    // Approximate p-value using empirical distribution
    // This is a simplified approximation
    if (test_stat < -4.0) return 0.001;  // Very strong rejection
    if (test_stat < -3.96) return 0.01;
    if (test_stat < -3.41) return 0.05;
    if (test_stat < -3.13) return 0.10;
    if (test_stat < -2.57) return 0.25;

    // For larger (less negative) values, p-value > 0.25
    // Use linear interpolation in the tail
    if (test_stat < -1.0) {
        // Linear interpolation between -2.57 (p=0.25) and -1.0 (p=0.90)
        double slope = (0.90 - 0.25) / (-1.0 + 2.57);
        return 0.25 + slope * (test_stat + 2.57);
    }

    return 0.99;  // Fail to reject null
}

// [[Rcpp::export]]
List cpp_adf_test(NumericVector x, int max_lag = -1, bool include_trend = true) {
    int n = x.size();

    // Default: use lag selection based on Schwert (1989)
    // lag = floor(12 * (n/100)^(1/4))
    if (max_lag < 0) {
        max_lag = std::floor(12.0 * std::pow(n / 100.0, 0.25));
        if (max_lag < 1) max_lag = 1;
    }

    // Ensure we have enough observations
    int min_obs = max_lag + 10;
    if (n < min_obs) {
        return List::create(
            Named("statistic") = NA_REAL,
            Named("p.value") = NA_REAL,
            Named("lags") = max_lag
        );
    }

    // Compute first differences
    int n_diff = n - 1;
    NumericVector diff(n_diff);
    for (int i = 0; i < n_diff; i++) {
        diff[i] = x[i + 1] - x[i];
    }

    // Set up regression with lags
    int start_idx = max_lag + 1;  // Start after we have enough lags
    int n_obs = n_diff - max_lag - 1;

    // Number of regressors: constant + trend (optional) + lagged level + lagged differences
    int n_reg = 2 + max_lag;  // constant + lagged level + max_lag lagged differences
    if (include_trend) n_reg++;

    NumericMatrix X(n_obs, n_reg);
    NumericVector y(n_obs);

    // Build regression matrices
    for (int i = 0; i < n_obs; i++) {
        int t = start_idx + i;  // Position in diff vector

        // Dependent variable: Δy_t
        y[i] = diff[t];

        int col = 0;

        // Constant
        X(i, col++) = 1.0;

        // Trend (if included)
        if (include_trend) {
            X(i, col++) = t + 1;  // Time index
        }

        // Lagged level: y_{t-1}
        X(i, col++) = x[t];  // This is y at position t in original series

        // Lagged differences: Δy_{t-1}, Δy_{t-2}, ..., Δy_{t-p}
        for (int lag = 1; lag <= max_lag; lag++) {
            X(i, col++) = diff[t - lag];
        }
    }

    // Run OLS regression
    List reg_result = ols_regression(X, y);
    NumericVector coefficients = reg_result["coefficients"];
    NumericVector std_errors = reg_result["std_errors"];

    // Check if regression succeeded
    if (NumericVector::is_na(coefficients[0])) {
        return List::create(
            Named("statistic") = NA_REAL,
            Named("p.value") = NA_REAL,
            Named("lags") = max_lag
        );
    }

    // Extract coefficient and std error for lagged level (γ coefficient)
    // Position: after constant and trend (if included)
    int gamma_idx = include_trend ? 2 : 1;
    double gamma = coefficients[gamma_idx];
    double se_gamma = std_errors[gamma_idx];

    // Compute t-statistic
    double t_stat = gamma / se_gamma;

    // Compute p-value using MacKinnon approximation
    double p_value = adf_pvalue_mackinnon(t_stat, n);

    return List::create(
        Named("statistic") = t_stat,
        Named("p.value") = p_value,
        Named("lags") = max_lag
    );
}
