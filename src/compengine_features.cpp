// CompEngine features from tsfeatures/hctsa
// These are features used in the CompEngine database
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <map>

using namespace Rcpp;

// Helper: First zero crossing of ACF (already exists, declared here for use)
int firstzero_ac_internal(const NumericVector& x);

// Helper: Compute ACF
NumericVector cpp_acf(NumericVector x, int max_lag, bool remove_mean);

//' @title Number of crossing points (median crossings)
//' @description Computes the number of times a time series crosses its median
//' @param x Numeric vector (time series)
//' @return Number of median crossings
//' @export
// [[Rcpp::export]]
int cpp_crossing_points(NumericVector x) {
  int n = x.size();
  if (n < 2) return 0;

  // Compute median
  NumericVector x_copy = clone(x);
  std::nth_element(x_copy.begin(), x_copy.begin() + n/2, x_copy.end());
  double midline = x_copy[n/2];
  if (n % 2 == 0) {
    std::nth_element(x_copy.begin(), x_copy.begin() + n/2 - 1, x_copy.end());
    midline = (midline + x_copy[n/2 - 1]) / 2.0;
  }

  // Count crossings
  int crossings = 0;
  for (int i = 0; i < n - 1; i++) {
    if (NumericVector::is_na(x[i]) || NumericVector::is_na(x[i+1])) continue;
    bool below_i = x[i] <= midline;
    bool below_ip1 = x[i+1] <= midline;
    if (below_i != below_ip1) {
      crossings++;
    }
  }

  return crossings;
}

//' @title Longest flat spot
//' @description Maximum run length within any of 10 equal-sized bins
//' @param x Numeric vector (time series)
//' @return Maximum run length in any bin
//' @export
// [[Rcpp::export]]
int cpp_flat_spots(NumericVector x) {
  int n = x.size();
  if (n < 2) return NA_INTEGER;

  // Remove NAs
  std::vector<double> x_clean;
  for (int i = 0; i < n; i++) {
    if (!NumericVector::is_na(x[i])) {
      x_clean.push_back(x[i]);
    }
  }

  if (x_clean.size() < 2) return NA_INTEGER;

  // Find min and max
  double x_min = *std::min_element(x_clean.begin(), x_clean.end());
  double x_max = *std::max_element(x_clean.begin(), x_clean.end());

  if (x_max == x_min) return static_cast<int>(x_clean.size()); // All same value

  // Create 10 bins
  double bin_width = (x_max - x_min) / 10.0;
  std::vector<int> bins(x_clean.size());

  for (size_t i = 0; i < x_clean.size(); i++) {
    int bin = static_cast<int>((x_clean[i] - x_min) / bin_width);
    if (bin >= 10) bin = 9; // Handle edge case where x[i] == x_max
    bins[i] = bin;
  }

  // Find maximum run length
  int max_run = 1;
  int current_run = 1;
  for (size_t i = 1; i < bins.size(); i++) {
    if (bins[i] == bins[i-1]) {
      current_run++;
      if (current_run > max_run) max_run = current_run;
    } else {
      current_run = 1;
    }
  }

  return max_run;
}

//' @title 2D embedding - proportion of points inside circle
//' @description Points inside given circular boundary in 2-d embedding space
//' @param y Numeric vector (time series)
//' @param boundary Circular boundary (1 or 2)
//' @return Proportion of points inside circle
//' @export
// [[Rcpp::export]]
double cpp_embed2_incircle(NumericVector y, int boundary) {
  int n = y.size();

  // Get first zero crossing of ACF for tau
  NumericVector acf_vals = cpp_acf(y, n-1, true);
  int tau = 0;
  for (int i = 1; i < acf_vals.size(); i++) {
    if (acf_vals[i] < 0) {
      tau = i;
      break;
    }
  }

  if (tau == 0 || tau >= n - 1) {
    return NA_REAL; // No zero crossing found or tau too large
  }

  // Create time-lagged vectors
  int N = n - tau;
  int count = 0;
  int valid = 0;

  for (int i = 0; i < N; i++) {
    double yt = y[i];
    double ytp = y[i + tau];

    if (!NumericVector::is_na(yt) && !NumericVector::is_na(ytp)) {
      valid++;
      if (yt*yt + ytp*ytp < boundary) {
        count++;
      }
    }
  }

  if (valid == 0) return NA_REAL;
  return static_cast<double>(count) / valid;
}

//' @title Binary motif entropy (3-letter words)
//' @description Entropy of 3-letter binary words (above/below mean)
//' @param y Numeric vector (time series)
//' @return Entropy of binary motifs
//' @export
// [[Rcpp::export]]
double cpp_motiftwo_entro3(NumericVector y) {
  int n = y.size();
  if (n < 5) {
    warning("Time series too short for motiftwo_entro3");
    return NA_REAL;
  }

  // Binarize: 1 if above mean, 0 if below
  double mean_y = mean(y);
  std::vector<int> yBin(n);
  for (int i = 0; i < n; i++) {
    yBin[i] = (y[i] > mean_y) ? 1 : 0;
  }

  // Count 3-letter patterns (000, 001, 010, 011, 100, 101, 110, 111)
  std::vector<int> counts(8, 0);
  for (int i = 0; i < n - 2; i++) {
    int pattern = yBin[i] * 4 + yBin[i+1] * 2 + yBin[i+2];
    counts[pattern]++;
  }

  // Compute probabilities and entropy
  int total = n - 2;
  double entropy = 0.0;
  for (int i = 0; i < 8; i++) {
    if (counts[i] > 0) {
      double p = static_cast<double>(counts[i]) / total;
      entropy -= p * std::log(p);
    }
  }

  return entropy;
}

//' @title Walker crossing proportion
//' @description Simulates a walker moving through time domain
//' @param y Numeric vector (time series)
//' @return Proportion of crossings between walker and time series
//' @export
// [[Rcpp::export]]
double cpp_walker_propcross(NumericVector y) {
  int n = y.size();
  if (n < 2) return NA_REAL;

  double p = 0.1; // Walker narrows gap by 10%
  std::vector<double> w(n);
  w[0] = 0.0; // Walker starts at zero

  // Walker follows: w[i] = w[i-1] + p * (y[i-1] - w[i-1])
  for (int i = 1; i < n; i++) {
    if (NumericVector::is_na(y[i-1])) {
      w[i] = w[i-1]; // Don't update if NA
    } else {
      w[i] = w[i-1] + p * (y[i-1] - w[i-1]);
    }
  }

  // Count crossings: (w[i] - y[i]) and (w[i+1] - y[i+1]) have opposite signs
  int crossings = 0;
  int valid = 0;
  for (int i = 0; i < n - 1; i++) {
    if (!NumericVector::is_na(y[i]) && !NumericVector::is_na(y[i+1])) {
      valid++;
      double diff_i = w[i] - y[i];
      double diff_ip1 = w[i+1] - y[i+1];
      if (diff_i * diff_ip1 < 0) {
        crossings++;
      }
    }
  }

  if (valid == 0) return NA_REAL;
  return static_cast<double>(crossings) / valid;
}
