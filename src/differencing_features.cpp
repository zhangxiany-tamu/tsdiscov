#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// Forward declarations from other files
NumericVector cpp_acf(NumericVector x, int max_lag, bool normalize);
NumericVector cpp_pacf(NumericVector x, int max_lag);

// [[Rcpp::export]]
NumericVector cpp_diff(NumericVector x, int lag = 1, int differences = 1) {
    // Compute differenced series: x[t] - x[t-lag]
    // Repeated 'differences' times

    int n = x.size();
    if (n < lag + 1) {
        return NumericVector(0);
    }

    NumericVector result = clone(x);

    for (int d = 0; d < differences; d++) {
        int len = result.size();
        if (len < lag + 1) break;

        NumericVector temp(len - lag);
        for (int i = 0; i < len - lag; i++) {
            temp[i] = result[i + lag] - result[i];
        }
        result = temp;
    }

    return result;
}

// [[Rcpp::export]]
double cpp_diff1_acf1(NumericVector x) {
    // First ACF of first differences

    if (x.size() < 3) return NA_REAL;

    NumericVector diff1 = cpp_diff(x, 1, 1);
    if (diff1.size() < 2) return NA_REAL;

    NumericVector acf = cpp_acf(diff1, 1, true);
    if (acf.size() < 1) return NA_REAL;

    return acf[0];
}

// [[Rcpp::export]]
double cpp_diff1_acf10(NumericVector x) {
    // Sum of squares of first 10 ACF lags of first differences

    if (x.size() < 12) return NA_REAL;

    NumericVector diff1 = cpp_diff(x, 1, 1);
    if (diff1.size() < 11) return NA_REAL;

    NumericVector acf = cpp_acf(diff1, 10, true);
    if (acf.size() < 10) return NA_REAL;

    double sum_sq = 0.0;
    for (int i = 0; i < 10; i++) {
        sum_sq += acf[i] * acf[i];
    }

    return sum_sq;
}

// [[Rcpp::export]]
double cpp_diff2_acf1(NumericVector x) {
    // First ACF of second differences

    if (x.size() < 4) return NA_REAL;

    NumericVector diff2 = cpp_diff(x, 1, 2);
    if (diff2.size() < 2) return NA_REAL;

    NumericVector acf = cpp_acf(diff2, 1, true);
    if (acf.size() < 1) return NA_REAL;

    return acf[0];
}

// [[Rcpp::export]]
double cpp_diff2_acf10(NumericVector x) {
    // Sum of squares of first 10 ACF lags of second differences

    if (x.size() < 13) return NA_REAL;

    NumericVector diff2 = cpp_diff(x, 1, 2);
    if (diff2.size() < 11) return NA_REAL;

    NumericVector acf = cpp_acf(diff2, 10, true);
    if (acf.size() < 10) return NA_REAL;

    double sum_sq = 0.0;
    for (int i = 0; i < 10; i++) {
        sum_sq += acf[i] * acf[i];
    }

    return sum_sq;
}

// [[Rcpp::export]]
double cpp_diff1x_pacf5(NumericVector x) {
    // Sum of squares of first 5 PACF lags of first differences

    if (x.size() < 7) return NA_REAL;

    NumericVector diff1 = cpp_diff(x, 1, 1);
    if (diff1.size() < 6) return NA_REAL;

    NumericVector pacf = cpp_pacf(diff1, 5);
    if (pacf.size() < 5) return NA_REAL;

    double sum_sq = 0.0;
    for (int i = 0; i < 5; i++) {
        sum_sq += pacf[i] * pacf[i];
    }

    return sum_sq;
}

// [[Rcpp::export]]
double cpp_diff2x_pacf5(NumericVector x) {
    // Sum of squares of first 5 PACF lags of second differences

    if (x.size() < 8) return NA_REAL;

    NumericVector diff2 = cpp_diff(x, 1, 2);
    if (diff2.size() < 6) return NA_REAL;

    NumericVector pacf = cpp_pacf(diff2, 5);
    if (pacf.size() < 5) return NA_REAL;

    double sum_sq = 0.0;
    for (int i = 0; i < 5; i++) {
        sum_sq += pacf[i] * pacf[i];
    }

    return sum_sq;
}

// [[Rcpp::export]]
NumericVector cpp_diff_acf_features(NumericVector x) {
    // Compute all differencing ACF features at once for efficiency

    NumericVector result = NumericVector::create(
        Named("diff1_acf1") = cpp_diff1_acf1(x),
        Named("diff1_acf10") = cpp_diff1_acf10(x),
        Named("diff2_acf1") = cpp_diff2_acf1(x),
        Named("diff2_acf10") = cpp_diff2_acf10(x),
        Named("diff1x_pacf5") = cpp_diff1x_pacf5(x),
        Named("diff2x_pacf5") = cpp_diff2x_pacf5(x)
    );

    return result;
}
