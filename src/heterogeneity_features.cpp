#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// Forward declaration
NumericVector cpp_acf(NumericVector x, int max_lag, bool normalize);

// [[Rcpp::export]]
double cpp_arch_acf(NumericVector x) {
    // Sum of squares of first 12 ACF of x²

    int n = x.size();
    if (n < 14) {
        return NA_REAL;
    }

    // Square the series
    NumericVector x_squared(n);
    for (int i = 0; i < n; i++) {
        x_squared[i] = x[i] * x[i];
    }

    // Compute ACF
    NumericVector acf_vals = cpp_acf(x_squared, 12, true);

    if (acf_vals.size() < 12) {
        return NA_REAL;
    }

    // Sum of squares
    double sum_sq = 0.0;
    for (int i = 0; i < 12; i++) {
        sum_sq += acf_vals[i] * acf_vals[i];
    }

    return sum_sq;
}
