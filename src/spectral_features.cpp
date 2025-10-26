#include <Rcpp.h>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// Helper function: Compute FFT using R's fft() function (much faster than naive DFT)
std::vector<std::complex<double>> compute_fft(const NumericVector& x) {
    int n = x.size();
    std::vector<std::complex<double>> result(n / 2 + 1);

    // Call R's fft() function for optimal performance
    Function fft("fft");
    ComplexVector r_fft = fft(x);

    // Extract first n/2 + 1 components (positive frequencies)
    for (int k = 0; k <= n / 2 && k < r_fft.size(); k++) {
        result[k] = std::complex<double>(r_fft[k].r, r_fft[k].i);
    }

    return result;
}

// Helper function: Compute power spectral density using Welch's method (simplified)
std::vector<double> welch_psd(const NumericVector& x) {
    int n = x.size();

    // For simplicity, use single segment (no windowing)
    // Full implementation would use overlapping windows
    std::vector<std::complex<double>> fft = compute_fft(x);
    std::vector<double> psd(fft.size());

    for (size_t i = 0; i < fft.size(); i++) {
        psd[i] = std::norm(fft[i]) / n;
    }

    return psd;
}

// [[Rcpp::export]]
double cpp_welch_power_low_freq(NumericVector x) {
    int n = x.size();
    if (n < 3) return NA_REAL;

    // Detrend by subtracting mean
    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= n;

    NumericVector detrended(n);
    for (int i = 0; i < n; i++) {
        detrended[i] = x[i] - mean;
    }

    // Compute power spectral density
    std::vector<double> psd = welch_psd(detrended);

    // Calculate power in lowest 20% of frequencies (indices 1 to 20% of spectrum)
    int low_freq_cutoff = std::max(1, static_cast<int>(psd.size() * 0.2));

    double total_power = 0.0;
    double low_freq_power = 0.0;

    for (size_t i = 1; i < psd.size(); i++) { // Skip DC component
        total_power += psd[i];
        if (static_cast<int>(i) <= low_freq_cutoff) {
            low_freq_power += psd[i];
        }
    }

    if (total_power == 0.0) return NA_REAL;
    return low_freq_power / total_power;
}

// [[Rcpp::export]]
double cpp_spectral_centroid(NumericVector x) {
    int n = x.size();
    if (n < 3) return NA_REAL;

    // Detrend by subtracting mean
    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= n;

    NumericVector detrended(n);
    for (int i = 0; i < n; i++) {
        detrended[i] = x[i] - mean;
    }

    // Compute power spectral density
    std::vector<double> psd = welch_psd(detrended);

    // Calculate spectral centroid (weighted mean of frequencies)
    double weighted_sum = 0.0;
    double total_power = 0.0;

    for (size_t i = 1; i < psd.size(); i++) { // Skip DC component
        double freq = static_cast<double>(i) / psd.size();
        weighted_sum += freq * psd[i];
        total_power += psd[i];
    }

    if (total_power == 0.0) return NA_REAL;
    return weighted_sum / total_power;
}

// [[Rcpp::export]]
double cpp_spectral_entropy(NumericVector x) {
    // Shannon entropy of power spectral density
    int n = x.size();
    if (n < 3) return NA_REAL;

    // Detrend by subtracting mean
    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= n;

    NumericVector detrended(n);
    for (int i = 0; i < n; i++) {
        detrended[i] = x[i] - mean;
    }

    // Compute power spectral density
    std::vector<double> psd = welch_psd(detrended);

    // Normalize PSD to get probability distribution
    double total_power = 0.0;
    for (size_t i = 1; i < psd.size(); i++) { // Skip DC
        total_power += psd[i];
    }

    if (total_power == 0.0) return NA_REAL;

    // Calculate Shannon entropy
    double entropy = 0.0;
    for (size_t i = 1; i < psd.size(); i++) {
        if (psd[i] > 0) {
            double p = psd[i] / total_power;
            entropy -= p * std::log(p);
        }
    }

    return entropy / std::log(2.0); // Normalize to bits
}

// [[Rcpp::export]]
NumericVector cpp_fft_coefficients(NumericVector x, int num_coef = 5) {
    // Extract FFT coefficients (magnitude, real, imaginary, angle)
    // for the first num_coef frequencies
    int n = x.size();
    if (n < 3 || num_coef < 1) {
        return NumericVector::create(Named("error") = NA_REAL);
    }

    // Limit num_coef to reasonable range
    num_coef = std::min(num_coef, n / 2);

    // Detrend by subtracting mean
    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= n;

    NumericVector detrended(n);
    for (int i = 0; i < n; i++) {
        detrended[i] = x[i] - mean;
    }

    // Compute FFT using R's fft()
    std::vector<std::complex<double>> fft = compute_fft(detrended);

    // Extract features for first num_coef frequencies (skip DC at k=0)
    NumericVector result(4 * num_coef);
    std::vector<std::string> names;

    for (int k = 1; k <= num_coef; k++) {
        double magnitude = std::abs(fft[k]);
        double real_part = fft[k].real();
        double imag_part = fft[k].imag();
        double angle = std::arg(fft[k]);

        int idx = (k - 1) * 4;
        result[idx] = magnitude;
        result[idx + 1] = real_part;
        result[idx + 2] = imag_part;
        result[idx + 3] = angle;

        names.push_back("fft_mag_" + std::to_string(k));
        names.push_back("fft_real_" + std::to_string(k));
        names.push_back("fft_imag_" + std::to_string(k));
        names.push_back("fft_angle_" + std::to_string(k));
    }

    result.names() = names;
    return result;
}

// [[Rcpp::export]]
NumericVector cpp_fft_aggregated(NumericVector x) {
    // Aggregate statistics of FFT magnitudes
    int n = x.size();
    if (n < 3) {
        return NumericVector::create(
            Named("fft_mean_mag") = NA_REAL,
            Named("fft_var_mag") = NA_REAL,
            Named("fft_max_mag") = NA_REAL
        );
    }

    // Detrend
    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= n;

    NumericVector detrended(n);
    for (int i = 0; i < n; i++) {
        detrended[i] = x[i] - mean;
    }

    // Compute FFT using R's fft()
    std::vector<std::complex<double>> fft = compute_fft(detrended);

    // Calculate statistics on magnitudes (skip DC at index 0)
    std::vector<double> magnitudes;
    for (size_t i = 1; i < fft.size(); i++) {
        magnitudes.push_back(std::abs(fft[i]));
    }

    if (magnitudes.empty()) {
        return NumericVector::create(
            Named("fft_mean_mag") = NA_REAL,
            Named("fft_var_mag") = NA_REAL,
            Named("fft_max_mag") = NA_REAL
        );
    }

    // Mean
    double sum = 0.0;
    for (double m : magnitudes) sum += m;
    double mean_mag = sum / magnitudes.size();

    // Variance
    double var_sum = 0.0;
    for (double m : magnitudes) {
        var_sum += (m - mean_mag) * (m - mean_mag);
    }
    double var_mag = var_sum / magnitudes.size();

    // Max
    double max_mag = *std::max_element(magnitudes.begin(), magnitudes.end());

    return NumericVector::create(
        Named("fft_mean_mag") = mean_mag,
        Named("fft_var_mag") = var_mag,
        Named("fft_max_mag") = max_mag
    );
}

// Helper function: Hanning window
std::vector<double> hanning_window(int n) {
    std::vector<double> window(n);
    for (int i = 0; i < n; i++) {
        window[i] = 0.5 * (1.0 - std::cos(2.0 * M_PI * i / (n - 1)));
    }
    return window;
}

// Welch's method for power spectral density estimation
// Returns PSD at n_freq evenly spaced frequencies
// [[Rcpp::export]]
NumericVector cpp_welch_psd(NumericVector x, int n_freq = 5, int window_size = -1, double overlap = 0.5) {
    int n = x.size();

    // Default result for edge cases
    NumericVector default_result(n_freq);
    for (int i = 0; i < n_freq; i++) {
        default_result[i] = NA_REAL;
    }

    // Need sufficient data
    if (n < 10) return default_result;

    // Remove NAs
    std::vector<double> x_clean;
    for (int i = 0; i < n; i++) {
        if (!NumericVector::is_na(x[i])) {
            x_clean.push_back(x[i]);
        }
    }

    if (x_clean.size() < 10) return default_result;
    n = x_clean.size();

    // Set default window size
    if (window_size < 0) {
        window_size = std::min(256, n / 2);
    }
    window_size = std::max(10, std::min(window_size, n));

    // Calculate step size based on overlap
    int step_size = std::max(1, static_cast<int>(window_size * (1.0 - overlap)));

    // Calculate number of segments
    int n_segments = 1 + (n - window_size) / step_size;
    if (n_segments < 1) {
        // Use single segment
        n_segments = 1;
        window_size = n;
        step_size = n;
    }

    // Create Hanning window
    std::vector<double> window = hanning_window(window_size);

    // Normalize window for energy preservation
    double window_norm = 0.0;
    for (int i = 0; i < window_size; i++) {
        window_norm += window[i] * window[i];
    }
    window_norm /= window_size;

    // Compute PSD for each segment and average
    std::vector<double> psd_avg(window_size / 2 + 1, 0.0);

    for (int seg = 0; seg < n_segments; seg++) {
        int start = seg * step_size;
        if (start + window_size > n) break;

        // Extract segment and apply window
        NumericVector segment(window_size);
        double seg_mean = 0.0;
        for (int i = 0; i < window_size; i++) {
            seg_mean += x_clean[start + i];
        }
        seg_mean /= window_size;

        for (int i = 0; i < window_size; i++) {
            segment[i] = (x_clean[start + i] - seg_mean) * window[i];
        }

        // Compute FFT of windowed segment
        std::vector<std::complex<double>> fft = compute_fft(segment);

        // Add power to average
        for (size_t i = 0; i < psd_avg.size(); i++) {
            psd_avg[i] += std::norm(fft[i]) / (window_size * window_norm);
        }
    }

    // Average over segments
    for (size_t i = 0; i < psd_avg.size(); i++) {
        psd_avg[i] /= n_segments;
    }

    // Extract power in n_freq bins
    NumericVector result(n_freq);
    int bins_per_freq = std::max(1, static_cast<int>(psd_avg.size() / n_freq));

    for (int f = 0; f < n_freq; f++) {
        double power = 0.0;
        int start_bin = f * bins_per_freq;
        int end_bin = (f == n_freq - 1) ? psd_avg.size() : (f + 1) * bins_per_freq;

        // Skip DC component (bin 0) in first frequency
        if (f == 0) start_bin = 1;

        for (int b = start_bin; b < end_bin && b < static_cast<int>(psd_avg.size()); b++) {
            power += psd_avg[b];
        }

        result[f] = power;
    }

    // Name the result
    CharacterVector names(n_freq);
    for (int i = 0; i < n_freq; i++) {
        names[i] = "welch_psd_bin" + std::to_string(i + 1);
    }
    result.names() = names;

    return result;
}
