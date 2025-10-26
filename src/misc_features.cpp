#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// Forward declaration from acf_features.cpp
NumericVector cpp_acf(NumericVector x, int max_lag, bool normalize);

// [[Rcpp::export]]
int cpp_firstzero_ac(NumericVector x, int max_lag = 50) {
    // Find first lag where ACF becomes negative
    // Returns the lag number (1-indexed), or 0 if no zero crossing found

    int n = x.size();
    if (n < 3) return 0;

    // Limit max_lag to reasonable range
    max_lag = std::min(max_lag, n - 1);

    // Compute ACF
    NumericVector acf_vals = cpp_acf(x, max_lag, true);

    // ACF returns values starting from lag 0, so skip lag 0
    // Find first negative value
    for (int i = 1; i < acf_vals.size(); i++) {
        if (acf_vals[i] < 0) {
            return i;  // Return lag number (1-indexed for lag 1, etc.)
        }
    }

    // No zero crossing found
    return 0;
}

// [[Rcpp::export]]
double cpp_zero_proportion(NumericVector x, double tol = 1e-8) {
    // Compute proportion of values close to zero

    int n = x.size();
    if (n == 0) return NA_REAL;

    int count_zeros = 0;
    int count_valid = 0;

    for (int i = 0; i < n; i++) {
        if (!NumericVector::is_na(x[i])) {
            count_valid++;
            if (std::abs(x[i]) < tol) {
                count_zeros++;
            }
        }
    }

    if (count_valid == 0) return NA_REAL;

    return static_cast<double>(count_zeros) / count_valid;
}

// [[Rcpp::export]]
double cpp_std1st_der(NumericVector x) {
    // Standard deviation of first differences

    int n = x.size();
    if (n < 2) return NA_REAL;

    // Compute first differences
    std::vector<double> diffs;
    for (int i = 1; i < n; i++) {
        if (!NumericVector::is_na(x[i]) && !NumericVector::is_na(x[i-1])) {
            diffs.push_back(x[i] - x[i-1]);
        }
    }

    if (diffs.size() < 2) return NA_REAL;

    // Compute mean of differences
    double mean_diff = 0.0;
    for (double d : diffs) {
        mean_diff += d;
    }
    mean_diff /= diffs.size();

    // Compute variance
    double var = 0.0;
    for (double d : diffs) {
        var += (d - mean_diff) * (d - mean_diff);
    }
    var /= (diffs.size() - 1);  // Sample variance

    return std::sqrt(var);
}
