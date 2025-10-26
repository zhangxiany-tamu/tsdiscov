#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
double cpp_absolute_sum_of_changes(NumericVector x) {
    // Sum of absolute differences
    int n = x.size();
    if (n < 2) return NA_REAL;

    double sum = 0.0;
    for (int i = 1; i < n; i++) {
        sum += std::abs(x[i] - x[i - 1]);
    }

    return sum;
}

// [[Rcpp::export]]
double cpp_range(NumericVector x) {
    // Range: max - min
    int n = x.size();
    if (n == 0) return NA_REAL;

    double min_val = *std::min_element(x.begin(), x.end());
    double max_val = *std::max_element(x.begin(), x.end());

    return max_val - min_val;
}

// [[Rcpp::export]]
double cpp_median_absolute_deviation(NumericVector x) {
    // Median absolute deviation from median
    int n = x.size();
    if (n == 0) return NA_REAL;

    // Calculate median
    NumericVector sorted = clone(x);
    std::sort(sorted.begin(), sorted.end());
    double median = (n % 2 == 0) ?
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0 :
        sorted[n / 2];

    // Calculate absolute deviations from median
    NumericVector abs_dev(n);
    for (int i = 0; i < n; i++) {
        abs_dev[i] = std::abs(x[i] - median);
    }

    // Return median of absolute deviations
    std::sort(abs_dev.begin(), abs_dev.end());
    double mad = (n % 2 == 0) ?
        (abs_dev[n / 2 - 1] + abs_dev[n / 2]) / 2.0 :
        abs_dev[n / 2];

    return mad;
}

// [[Rcpp::export]]
double cpp_coefficient_of_variation(NumericVector x) {
    // Coefficient of variation: std / mean
    int n = x.size();
    if (n < 2) return NA_REAL;

    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= n;

    if (mean == 0.0) return NA_REAL;

    double variance = 0.0;
    for (int i = 0; i < n; i++) {
        double diff = x[i] - mean;
        variance += diff * diff;
    }
    double std = std::sqrt(variance / (n - 1));

    return std / std::abs(mean);
}

// [[Rcpp::export]]
double cpp_benford_correlation(NumericVector x) {
    // Correlation with Benford's Law distribution
    int n = x.size();
    if (n < 10) return NA_REAL;

    // Expected Benford proportions for first digits 1-9
    std::vector<double> benford_prob = {
        0.30103, 0.17609, 0.12494, 0.09691, 0.07918,
        0.06695, 0.05799, 0.05115, 0.04576
    };

    // Count first digits in absolute values
    std::vector<int> digit_counts(9, 0);
    int valid_count = 0;

    for (int i = 0; i < n; i++) {
        double abs_val = std::abs(x[i]);
        if (abs_val > 1e-10) { // Avoid issues with very small values
            // Extract first digit
            while (abs_val >= 10.0) abs_val /= 10.0;
            while (abs_val < 1.0) abs_val *= 10.0;

            int first_digit = static_cast<int>(abs_val);
            if (first_digit >= 1 && first_digit <= 9) {
                digit_counts[first_digit - 1]++;
                valid_count++;
            }
        }
    }

    if (valid_count < 10) return NA_REAL;

    // Calculate observed proportions
    std::vector<double> observed_prob(9);
    for (int i = 0; i < 9; i++) {
        observed_prob[i] = static_cast<double>(digit_counts[i]) / valid_count;
    }

    // Calculate Pearson correlation
    double mean_benford = 0.0, mean_observed = 0.0;
    for (int i = 0; i < 9; i++) {
        mean_benford += benford_prob[i];
        mean_observed += observed_prob[i];
    }
    mean_benford /= 9.0;
    mean_observed /= 9.0;

    double numerator = 0.0, var_benford = 0.0, var_observed = 0.0;
    for (int i = 0; i < 9; i++) {
        double diff_b = benford_prob[i] - mean_benford;
        double diff_o = observed_prob[i] - mean_observed;
        numerator += diff_b * diff_o;
        var_benford += diff_b * diff_b;
        var_observed += diff_o * diff_o;
    }

    if (var_benford == 0.0 || var_observed == 0.0) return NA_REAL;

    return numerator / std::sqrt(var_benford * var_observed);
}
