#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
double cpp_mean(NumericVector x) {
    int n = x.size();
    if (n == 0) return NA_REAL;
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        if (!R_IsNA(x[i])) sum += x[i];
    }
    return sum / n;
}

// [[Rcpp::export]]
double cpp_std(NumericVector x) {
    int n = x.size();
    if (n < 2) return NA_REAL;

    double mean = cpp_mean(x);
    double sum_sq = 0.0;

    for (int i = 0; i < n; i++) {
        if (!R_IsNA(x[i])) {
            double diff = x[i] - mean;
            sum_sq += diff * diff;
        }
    }
    return std::sqrt(sum_sq / (n - 1));
}

// [[Rcpp::export]]
double cpp_variance(NumericVector x) {
    double s = cpp_std(x);
    return s * s;
}

// [[Rcpp::export]]
double cpp_skewness(NumericVector x) {
    int n = x.size();
    if (n < 3) return NA_REAL;

    double mean = cpp_mean(x);
    double std = cpp_std(x);
    if (std == 0) return NA_REAL;

    double sum_cube = 0.0;
    for (int i = 0; i < n; i++) {
        if (!R_IsNA(x[i])) {
            double z = (x[i] - mean) / std;
            sum_cube += z * z * z;
        }
    }
    return sum_cube / n;
}

// [[Rcpp::export]]
double cpp_kurtosis(NumericVector x) {
    int n = x.size();
    if (n < 4) return NA_REAL;

    double mean = cpp_mean(x);
    double std = cpp_std(x);
    if (std == 0) return NA_REAL;

    double sum_fourth = 0.0;
    for (int i = 0; i < n; i++) {
        if (!R_IsNA(x[i])) {
            double z = (x[i] - mean) / std;
            sum_fourth += z * z * z * z;
        }
    }
    return sum_fourth / n - 3.0;  // Excess kurtosis
}

// [[Rcpp::export]]
NumericVector cpp_quantile(NumericVector x, NumericVector probs) {
    NumericVector sorted = clone(x);
    std::sort(sorted.begin(), sorted.end());

    int n = sorted.size();
    int np = probs.size();
    NumericVector result(np);

    for (int i = 0; i < np; i++) {
        double h = (n - 1) * probs[i];
        int h_floor = static_cast<int>(std::floor(h));
        double w = h - h_floor;

        if (h_floor >= n - 1) {
            result[i] = sorted[n - 1];
        } else {
            result[i] = (1 - w) * sorted[h_floor] + w * sorted[h_floor + 1];
        }
    }
    return result;
}
