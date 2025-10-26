#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

//' Compute spectral shape descriptors
//'
//' @param x Numeric vector
//' @return Named list with spectral shape features
//' @export
// [[Rcpp::export]]
List cpp_spectral_shape(NumericVector x) {
  int n = x.size();

  // Default return
  List default_result = List::create(
    Named("spectral_flatness") = NA_REAL,
    Named("spectral_slope") = NA_REAL,
    Named("spectral_rolloff_95") = NA_REAL,
    Named("spectral_bandwidth") = NA_REAL
  );

  if (n < 20) {
    return default_result;
  }

  // Compute FFT magnitude (power spectrum)
  int nfft = n;
  int nhalf = nfft / 2 + 1;

  std::vector<std::complex<double>> fft_result(nfft);

  // Simple DFT (for positive frequencies only, we'll use magnitude)
  std::vector<double> psd(nhalf);

  for (int k = 0; k < nhalf; k++) {
    double real_part = 0.0;
    double imag_part = 0.0;

    for (int j = 0; j < n; j++) {
      double angle = -2.0 * M_PI * k * j / n;
      real_part += x[j] * std::cos(angle);
      imag_part += x[j] * std::sin(angle);
    }

    // Power spectral density (magnitude squared)
    psd[k] = std::sqrt(real_part * real_part + imag_part * imag_part) / n;
  }

  // Avoid division by zero
  double min_psd = 1e-10;
  for (int i = 0; i < nhalf; i++) {
    if (psd[i] < min_psd) psd[i] = min_psd;
  }

  // 1. Spectral flatness (geometric mean / arithmetic mean)
  double geo_mean_log = 0.0;
  double arith_mean = 0.0;

  for (int i = 0; i < nhalf; i++) {
    geo_mean_log += std::log(psd[i]);
    arith_mean += psd[i];
  }

  geo_mean_log /= nhalf;
  arith_mean /= nhalf;

  double spectral_flatness = std::exp(geo_mean_log) / arith_mean;

  // 2. Spectral slope (linear regression on log spectrum)
  std::vector<double> freqs(nhalf);
  for (int i = 0; i < nhalf; i++) {
    freqs[i] = (double)i / nhalf;
  }

  double mean_freq = 0.0;
  double mean_log_psd = 0.0;
  for (int i = 0; i < nhalf; i++) {
    mean_freq += freqs[i];
    mean_log_psd += std::log(psd[i]);
  }
  mean_freq /= nhalf;
  mean_log_psd /= nhalf;

  double numerator = 0.0;
  double denominator = 0.0;
  for (int i = 0; i < nhalf; i++) {
    numerator += (freqs[i] - mean_freq) * (std::log(psd[i]) - mean_log_psd);
    denominator += (freqs[i] - mean_freq) * (freqs[i] - mean_freq);
  }

  double spectral_slope = (denominator > 1e-10) ? (numerator / denominator) : 0.0;

  // 3. Spectral rolloff (95% energy threshold)
  double total_energy = 0.0;
  for (int i = 0; i < nhalf; i++) {
    total_energy += psd[i] * psd[i];
  }

  double threshold = 0.95 * total_energy;
  double cumulative_energy = 0.0;
  int rolloff_idx = nhalf - 1;

  for (int i = 0; i < nhalf; i++) {
    cumulative_energy += psd[i] * psd[i];
    if (cumulative_energy >= threshold) {
      rolloff_idx = i;
      break;
    }
  }

  double spectral_rolloff_95 = (double)rolloff_idx / nhalf;

  // 4. Spectral bandwidth (standard deviation weighted by magnitude)
  double weighted_mean = 0.0;
  double total_weight = 0.0;

  for (int i = 0; i < nhalf; i++) {
    weighted_mean += freqs[i] * psd[i];
    total_weight += psd[i];
  }

  weighted_mean /= (total_weight > 1e-10 ? total_weight : 1.0);

  double variance = 0.0;
  for (int i = 0; i < nhalf; i++) {
    double diff = freqs[i] - weighted_mean;
    variance += diff * diff * psd[i];
  }
  variance /= (total_weight > 1e-10 ? total_weight : 1.0);

  double spectral_bandwidth = std::sqrt(variance);

  // Return results
  return List::create(
    Named("spectral_flatness") = spectral_flatness,
    Named("spectral_slope") = spectral_slope,
    Named("spectral_rolloff_95") = spectral_rolloff_95,
    Named("spectral_bandwidth") = spectral_bandwidth
  );
}
