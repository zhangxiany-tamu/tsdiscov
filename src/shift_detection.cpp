#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// Helper: Compute rolling mean efficiently using sliding window
std::vector<double> rolling_mean(const NumericVector& x, int width) {
    int n = x.size();
    std::vector<double> result;

    if (n < width) {
        return result;
    }

    // Compute first window
    double sum = 0.0;
    for (int i = 0; i < width; i++) {
        sum += x[i];
    }
    result.push_back(sum / width);

    // Slide window
    for (int i = width; i < n; i++) {
        sum = sum - x[i - width] + x[i];
        result.push_back(sum / width);
    }

    return result;
}

// Helper: Compute rolling variance efficiently
std::vector<double> rolling_var(const NumericVector& x, int width) {
    int n = x.size();
    std::vector<double> result;

    if (n < width) {
        return result;
    }

    // For each window, compute variance
    for (int i = 0; i <= n - width; i++) {
        double sum = 0.0;
        double sum_sq = 0.0;

        for (int j = i; j < i + width; j++) {
            sum += x[j];
            sum_sq += x[j] * x[j];
        }

        double mean = sum / width;
        double variance = (sum_sq / width) - (mean * mean);
        result.push_back(variance);
    }

    return result;
}

// Helper: Silverman's rule of thumb for bandwidth selection
double bandwidth_nrd0(const NumericVector& x) {
    int n = x.size();
    if (n < 2) return 1.0;

    // Compute standard deviation
    double mean = 0.0;
    for (int i = 0; i < n; i++) {
        mean += x[i];
    }
    mean /= n;

    double var = 0.0;
    for (int i = 0; i < n; i++) {
        var += (x[i] - mean) * (x[i] - mean);
    }
    var /= (n - 1);
    double sd = std::sqrt(var);

    // Compute IQR
    NumericVector sorted = clone(x);
    std::sort(sorted.begin(), sorted.end());
    double q25 = sorted[n / 4];
    double q75 = sorted[3 * n / 4];
    double iqr = q75 - q25;

    // Silverman's rule: bw = 0.9 * min(sd, IQR/1.34) * n^(-1/5)
    double lo = std::min(sd, iqr / 1.34);
    if (lo == 0) lo = sd;
    if (lo == 0) lo = 1.0;

    return 0.9 * lo * std::pow(n, -0.2);
}

// Helper: Gaussian kernel density
double dnorm_kernel(double x, double mean, double sd) {
    const double pi = 3.14159265358979323846;
    double z = (x - mean) / sd;
    return std::exp(-0.5 * z * z) / (sd * std::sqrt(2.0 * pi));
}

// [[Rcpp::export]]
NumericVector cpp_shift_detection(NumericVector x, int width = 10) {
    // Compute all 6 shift detection features
    // Returns: max_level_shift, time_level_shift, max_var_shift,
    //          time_var_shift, max_kl_shift, time_kl_shift

    int n = x.size();

    // Handle short series
    if (n < 2 * width) {
        return NumericVector::create(
            Named("max_level_shift") = NA_REAL,
            Named("time_level_shift") = NA_REAL,
            Named("max_var_shift") = NA_REAL,
            Named("time_var_shift") = NA_REAL,
            Named("max_kl_shift") = NA_REAL,
            Named("time_kl_shift") = NA_REAL
        );
    }

    // 1. Level shift (mean shift)
    std::vector<double> rmean = rolling_mean(x, width);
    double max_level = 0.0;
    int time_level = NA_INTEGER;

    if (rmean.size() > static_cast<size_t>(width)) {
        for (size_t i = width; i < rmean.size(); i++) {
            double shift = std::abs(rmean[i] - rmean[i - width]);
            if (shift > max_level) {
                max_level = shift;
                time_level = i + width - 1;  // Adjust for R's 1-indexing
            }
        }
    }

    // 2. Variance shift
    std::vector<double> rvar = rolling_var(x, width);
    double max_var = 0.0;
    int time_var = NA_INTEGER;

    if (rvar.size() > static_cast<size_t>(width)) {
        for (size_t i = width; i < rvar.size(); i++) {
            double shift = std::abs(rvar[i] - rvar[i - width]);
            if (shift > max_var) {
                max_var = shift;
                time_var = i + width - 1;
            }
        }
    }

    // 3. KL divergence shift
    double max_kl = 0.0;
    int time_kl = NA_INTEGER;

    // Remove NAs for bandwidth calculation
    NumericVector x_clean;
    for (int i = 0; i < n; i++) {
        if (!NumericVector::is_na(x[i])) {
            x_clean.push_back(x[i]);
        }
    }

    if (static_cast<int>(x_clean.size()) >= 2 * width) {
        double bw = bandwidth_nrd0(x_clean);

        // Create grid for density estimation
        double xmin = *std::min_element(x_clean.begin(), x_clean.end());
        double xmax = *std::max_element(x_clean.begin(), x_clean.end());
        int gw = 100;  // grid width
        std::vector<double> xgrid(gw);
        double grid_step = (xmax - xmin) / (gw - 1);
        for (int i = 0; i < gw; i++) {
            xgrid[i] = xmin + i * grid_step;
        }

        // Compute density matrix (each row is density at time i)
        std::vector<std::vector<double>> dens_mat(n, std::vector<double>(gw));
        double min_density = dnorm_kernel(38, 0, 1);  // Very small probability

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < gw; j++) {
                double d = dnorm_kernel(xgrid[j], x[i], bw);
                dens_mat[i][j] = std::max(d, min_density);
            }
        }

        // Compute rolling mean of density (by column)
        std::vector<std::vector<double>> rmean_dens(n, std::vector<double>(gw, 0.0));
        for (int col = 0; col < gw; col++) {
            double sum = 0.0;
            for (int i = 0; i < width; i++) {
                sum += dens_mat[i][col];
            }
            rmean_dens[width - 1][col] = sum / width;

            for (int i = width; i < n; i++) {
                sum = sum - dens_mat[i - width][col] + dens_mat[i][col];
                rmean_dens[i][col] = sum / width;
            }
        }

        // Compute KL divergence between consecutive windows
        int max_idx = n - width;
        for (int i = 0; i < max_idx; i++) {
            int lo = i;
            int hi = i + width;

            if (hi >= n) break;

            // KL divergence: sum(p * log(p/q) * grid_step)
            double kl = 0.0;
            for (int j = 0; j < gw; j++) {
                double p = rmean_dens[lo][j];
                double q = rmean_dens[hi][j];
                if (p > min_density && q > min_density) {
                    kl += p * std::log(p / q) * grid_step;
                }
            }

            if (kl > max_kl) {
                max_kl = kl;
                time_kl = hi;
            }
        }
    }

    return NumericVector::create(
        Named("max_level_shift") = max_level,
        Named("time_level_shift") = time_level,
        Named("max_var_shift") = max_var,
        Named("time_var_shift") = time_var,
        Named("max_kl_shift") = max_kl,
        Named("time_kl_shift") = time_kl
    );
}
