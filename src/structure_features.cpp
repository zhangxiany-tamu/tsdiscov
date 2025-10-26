#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
int cpp_number_peaks(NumericVector x, int support = 3) {
    // Count number of peaks (local maxima)
    int n = x.size();
    if (n < 3) return 0;

    int peak_count = 0;

    for (int i = support; i < n - support; i++) {
        bool is_peak = true;

        // Check if x[i] is greater than all points in support window
        for (int j = i - support; j <= i + support; j++) {
            if (j != i && x[j] >= x[i]) {
                is_peak = false;
                break;
            }
        }

        if (is_peak) {
            peak_count++;
        }
    }

    return peak_count;
}

// [[Rcpp::export]]
int cpp_longest_strike_above_mean(NumericVector x) {
    // Longest consecutive run of values above mean
    int n = x.size();
    if (n == 0) return 0;

    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= n;

    int max_strike = 0;
    int current_strike = 0;

    for (int i = 0; i < n; i++) {
        if (x[i] > mean) {
            current_strike++;
            if (current_strike > max_strike) {
                max_strike = current_strike;
            }
        } else {
            current_strike = 0;
        }
    }

    return max_strike;
}

// [[Rcpp::export]]
int cpp_longest_strike_below_mean(NumericVector x) {
    // Longest consecutive run of values below mean
    int n = x.size();
    if (n == 0) return 0;

    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= n;

    int max_strike = 0;
    int current_strike = 0;

    for (int i = 0; i < n; i++) {
        if (x[i] < mean) {
            current_strike++;
            if (current_strike > max_strike) {
                max_strike = current_strike;
            }
        } else {
            current_strike = 0;
        }
    }

    return max_strike;
}

// [[Rcpp::export]]
int cpp_count_above_mean(NumericVector x) {
    // Count values above mean
    int n = x.size();
    if (n == 0) return 0;

    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= n;

    int count = 0;
    for (int i = 0; i < n; i++) {
        if (x[i] > mean) count++;
    }

    return count;
}

// [[Rcpp::export]]
int cpp_count_below_mean(NumericVector x) {
    // Count values below mean
    int n = x.size();
    if (n == 0) return 0;

    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= n;

    int count = 0;
    for (int i = 0; i < n; i++) {
        if (x[i] < mean) count++;
    }

    return count;
}

// [[Rcpp::export]]
double cpp_ratio_beyond_r_sigma(NumericVector x, double r = 1.0) {
    // Ratio of values beyond r standard deviations from mean
    int n = x.size();
    if (n < 2) return NA_REAL;

    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= n;

    double variance = 0.0;
    for (int i = 0; i < n; i++) {
        double diff = x[i] - mean;
        variance += diff * diff;
    }
    double std = std::sqrt(variance / (n - 1));

    if (std == 0.0) return 0.0;

    double threshold = r * std;
    int count = 0;

    for (int i = 0; i < n; i++) {
        if (std::abs(x[i] - mean) > threshold) {
            count++;
        }
    }

    return static_cast<double>(count) / n;
}

// [[Rcpp::export]]
int cpp_number_crossings(NumericVector x, double m = 0.0) {
    // Number of times series crosses value m
    int n = x.size();
    if (n < 2) return 0;

    int crossings = 0;
    bool above = (x[0] > m);

    for (int i = 1; i < n; i++) {
        bool current_above = (x[i] > m);
        if (current_above != above) {
            crossings++;
            above = current_above;
        }
    }

    return crossings;
}
