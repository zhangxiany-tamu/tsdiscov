#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Compute pairwise cross-correlation function (CCF) features - Rcpp version
//'
//' Efficiently computes CCF for all pairs of time series
//'
//' @param X Matrix (N x T) where rows are series, columns are timepoints
//' @param max_lag Maximum lag for CCF
//' @return List with ccf_max_vals, ccf_lag_vals, ccf_zero_vals for all pairs
//' @keywords internal
// [[Rcpp::export]]
List cpp_pairwise_ccf(arma::mat X, int max_lag) {
    int N = X.n_rows;
    int T = X.n_cols;
    int n_pairs = N * (N - 1) / 2;

    NumericVector ccf_max_vals(n_pairs);
    IntegerVector ccf_lag_vals(n_pairs);
    NumericVector ccf_zero_vals(n_pairs);

    int pair_idx = 0;

    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            // Extract series
            arma::rowvec x = X.row(i);
            arma::rowvec y = X.row(j);

            // Demean
            double mean_x = arma::mean(x);
            double mean_y = arma::mean(y);
            arma::rowvec x_dm = x - mean_x;
            arma::rowvec y_dm = y - mean_y;

            // Standard deviations
            double sd_x = arma::stddev(x);
            double sd_y = arma::stddev(y);

            if (sd_x < 1e-10 || sd_y < 1e-10) {
                ccf_max_vals[pair_idx] = NA_REAL;
                ccf_lag_vals[pair_idx] = 0;
                ccf_zero_vals[pair_idx] = NA_REAL;
                pair_idx++;
                continue;
            }

            // Compute CCF for lags from -max_lag to +max_lag
            int n_lags = 2 * max_lag + 1;
            arma::vec ccf_vals(n_lags);

            for (int lag = -max_lag; lag <= max_lag; lag++) {
                double sum = 0.0;
                int count = 0;

                // Compute correlation at this lag
                for (int t = 0; t < T; t++) {
                    int t_lag = t + lag;
                    if (t_lag >= 0 && t_lag < T) {
                        sum += x_dm[t] * y_dm[t_lag];
                        count++;
                    }
                }

                ccf_vals[lag + max_lag] = (count > 0) ? sum / (count * sd_x * sd_y) : 0.0;
            }

            // Find maximum absolute CCF
            arma::vec abs_ccf = arma::abs(ccf_vals);
            arma::uword max_idx = abs_ccf.index_max();
            ccf_max_vals[pair_idx] = abs_ccf[max_idx];
            ccf_lag_vals[pair_idx] = static_cast<int>(max_idx) - max_lag;

            // CCF at zero lag
            ccf_zero_vals[pair_idx] = std::abs(ccf_vals[max_lag]);

            pair_idx++;
        }
    }

    return List::create(
        Named("ccf_max_vals") = ccf_max_vals,
        Named("ccf_lag_vals") = ccf_lag_vals,
        Named("ccf_zero_vals") = ccf_zero_vals
    );
}

//' Compute mutual information between two vectors - Rcpp version
//'
//' Histogram-based estimator of mutual information
//'
//' @param x Numeric vector
//' @param y Numeric vector
//' @param bins Number of bins for histogram
//' @return Mutual information estimate
//' @keywords internal
// [[Rcpp::export]]
double cpp_mutual_information(NumericVector x, NumericVector y, int bins = 10) {
    int n = x.size();
    if (n < 10) return NA_REAL;

    // Find ranges
    double x_min = *std::min_element(x.begin(), x.end());
    double x_max = *std::max_element(x.begin(), x.end());
    double y_min = *std::min_element(y.begin(), y.end());
    double y_max = *std::max_element(y.begin(), y.end());

    if (x_max - x_min < 1e-10 || y_max - y_min < 1e-10) {
        return 0.0;
    }

    double x_width = (x_max - x_min) / bins;
    double y_width = (y_max - y_min) / bins;

    // Create joint histogram
    arma::mat joint_hist(bins, bins, arma::fill::zeros);

    for (int i = 0; i < n; i++) {
        int x_bin = std::min(static_cast<int>((x[i] - x_min) / x_width), bins - 1);
        int y_bin = std::min(static_cast<int>((y[i] - y_min) / y_width), bins - 1);
        joint_hist(x_bin, y_bin) += 1.0;
    }

    // Normalize to get probabilities
    arma::mat p_xy = joint_hist / n;

    // Marginal probabilities
    arma::vec p_x = arma::sum(p_xy, 1);
    arma::rowvec p_y = arma::sum(p_xy, 0);

    // Compute mutual information
    double mi = 0.0;
    for (int i = 0; i < bins; i++) {
        for (int j = 0; j < bins; j++) {
            if (p_xy(i, j) > 1e-10 && p_x[i] > 1e-10 && p_y[j] > 1e-10) {
                mi += p_xy(i, j) * std::log(p_xy(i, j) / (p_x[i] * p_y[j]));
            }
        }
    }

    return mi;
}

//' Compute pairwise mutual information - Rcpp version
//'
//' Efficiently computes MI for all pairs of time series
//'
//' @param X Matrix (N x T) where rows are series, columns are timepoints
//' @param bins Number of bins for histogram
//' @return Vector of mutual information values for all pairs
//' @keywords internal
// [[Rcpp::export]]
NumericVector cpp_pairwise_mi(arma::mat X, int bins = 10) {
    int N = X.n_rows;
    int n_pairs = N * (N - 1) / 2;

    NumericVector mi_vals(n_pairs);
    int pair_idx = 0;

    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            NumericVector x = as<NumericVector>(wrap(X.row(i)));
            NumericVector y = as<NumericVector>(wrap(X.row(j)));
            mi_vals[pair_idx] = cpp_mutual_information(x, y, bins);
            pair_idx++;
        }
    }

    return mi_vals;
}

//' Compute Shannon entropy from binned data - multivariate version
//'
//' @param x Numeric vector
//' @param bins Number of bins
//' @return Shannon entropy
//' @keywords internal
// [[Rcpp::export]]
double cpp_mv_shannon_entropy(NumericVector x, int bins = 10) {
    int n = x.size();
    if (n < 10) return NA_REAL;

    double x_min = *std::min_element(x.begin(), x.end());
    double x_max = *std::max_element(x.begin(), x.end());

    if (x_max - x_min < 1e-10) return 0.0;

    double width = (x_max - x_min) / bins;

    // Create histogram
    std::vector<int> hist(bins, 0);
    for (int i = 0; i < n; i++) {
        int bin = std::min(static_cast<int>((x[i] - x_min) / width), bins - 1);
        hist[bin]++;
    }

    // Compute entropy
    double entropy = 0.0;
    for (int i = 0; i < bins; i++) {
        if (hist[i] > 0) {
            double p = static_cast<double>(hist[i]) / n;
            entropy -= p * std::log(p);
        }
    }

    return entropy;
}

//' Compute cross-spectral density for a pair of series - Rcpp version
//'
//' Uses FFT to compute cross-spectral density
//'
//' @param x First time series
//' @param y Second time series
//' @return List with coherence, phase, and CSD magnitude
//' @keywords internal
// [[Rcpp::export]]
List cpp_cross_spectral_density(NumericVector x, NumericVector y) {
    int n = x.size();
    if (n < 4) {
        return List::create(
            Named("coherence") = NumericVector(0),
            Named("phase") = NumericVector(0),
            Named("csd_mag") = 0.0
        );
    }

    // Demean
    double mean_x = std::accumulate(x.begin(), x.end(), 0.0) / n;
    double mean_y = std::accumulate(y.begin(), y.end(), 0.0) / n;

    // Prepare for FFT (using R's FFT from stats)
    // Note: This is a placeholder - in practice we'd use FFTW or similar
    // For now, we'll return to R for FFT computation and just provide helper functions

    return List::create(
        Named("note") = "FFT computation delegated to R's stats::fft"
    );
}

//' Compute per-series statistics for diversity features - Rcpp version
//'
//' Efficiently computes univariate statistics for each series
//'
//' @param X Matrix (N x T) where rows are series
//' @return List with means, variances, ranges, skewness, kurtosis for each series
//' @keywords internal
// [[Rcpp::export]]
List cpp_series_statistics(arma::mat X) {
    int N = X.n_rows;

    NumericVector means(N);
    NumericVector variances(N);
    NumericVector ranges(N);
    NumericVector skewness(N);
    NumericVector kurtosis(N);

    for (int i = 0; i < N; i++) {
        arma::rowvec series = X.row(i);

        // Mean
        double m = arma::mean(series);
        means[i] = m;

        // Variance
        double v = arma::var(series);
        variances[i] = v;

        // Range
        ranges[i] = arma::max(series) - arma::min(series);

        // Standardize for higher moments
        if (v > 1e-10) {
            arma::rowvec z = (series - m) / std::sqrt(v);

            // Skewness
            skewness[i] = arma::mean(arma::pow(z, 3));

            // Excess kurtosis
            kurtosis[i] = arma::mean(arma::pow(z, 4)) - 3.0;
        } else {
            skewness[i] = 0.0;
            kurtosis[i] = 0.0;
        }
    }

    return List::create(
        Named("means") = means,
        Named("variances") = variances,
        Named("ranges") = ranges,
        Named("skewness") = skewness,
        Named("kurtosis") = kurtosis
    );
}
