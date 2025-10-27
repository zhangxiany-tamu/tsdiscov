#' Cross-spectral and coherence features for multivariate time series
#'
#' Extracts features characterizing frequency-domain relationships between series.
#' Uses magnitude-squared coherence and cross-spectral density to measure
#' frequency-specific synchronization.
#'
#' Note: This function is exported for direct use but is not currently included
#' in the default "all" feature set or accessible via ts_features_multivariate().
#' Call this function directly: ts_mv_spectral(X).
#'
#' @param X Matrix (N x T) where N is number of series, T is timepoints
#' @return Named list of 10 cross-spectral features
#' @export
#' @examples
#' \dontrun{
#' # Create signals with coherent low-frequency component
#' t <- seq(0, 10, length.out = 256)
#' common <- sin(2*pi*t) + sin(4*pi*t)
#' X <- matrix(0, nrow = 3, ncol = 256)
#' X[1,] <- common + rnorm(256, sd=0.1)
#' X[2,] <- common + rnorm(256, sd=0.1)
#' X[3,] <- rnorm(256)  # Independent
#' spectral_features <- ts_mv_spectral(X)
#' }
ts_mv_spectral <- function(X) {
  N <- nrow(X)
  T_len <- ncol(X)

  if (N < 2) {
    # Single series - no cross-spectral analysis
    return(list(
      coherence_mean = NA, coherence_median = NA, coherence_max = NA,
      coherence_std = NA, coherence_high_frac = NA,
      phase_lag_mean = NA, phase_lag_std = NA,
      spectral_corr_mean = NA, spectral_corr_max = NA,
      csd_magnitude_mean = NA
    ))
  }

  # Compute power spectral density for each series using Welch's method
  # Use simple periodogram for simplicity (can enhance with Welch later)
  psds <- list()
  for (i in 1:N) {
    spec_result <- tryCatch({
      spectrum(X[i,], plot = FALSE, method = "pgram")
    }, error = function(e) NULL)

    if (!is.null(spec_result)) {
      psds[[i]] <- list(
        freq = spec_result$freq,
        spec = spec_result$spec
      )
    } else {
      psds[[i]] <- NULL
    }
  }

  # Check if all PSDs computed successfully
  if (any(sapply(psds, is.null))) {
    return(list(
      coherence_mean = NA, coherence_median = NA, coherence_max = NA,
      coherence_std = NA, coherence_high_frac = NA,
      phase_lag_mean = NA, phase_lag_std = NA,
      spectral_corr_mean = NA, spectral_corr_max = NA,
      csd_magnitude_mean = NA
    ))
  }

  # Compute coherence and phase for all pairs
  n_pairs <- N * (N - 1) / 2
  coherence_vals <- list()
  phase_lags <- numeric(n_pairs)
  csd_mags <- numeric(n_pairs)
  spectral_corrs <- numeric(n_pairs)

  pair_idx <- 1
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      # Compute cross-spectrum using FFT
      fft_i <- fft(X[i,])
      fft_j <- fft(X[j,])

      # Cross-spectral density
      csd <- fft_i * Conj(fft_j)

      # Use only positive frequencies (first half + Nyquist)
      n_freq <- floor(T_len / 2) + 1
      csd_pos <- csd[1:n_freq]

      # Auto-spectra
      psd_i <- Mod(fft_i[1:n_freq])^2
      psd_j <- Mod(fft_j[1:n_freq])^2

      # Magnitude-squared coherence: |Pxy|^2 / (Pxx * Pyy)
      coh <- (Mod(csd_pos)^2) / (psd_i * psd_j + 1e-10)
      coh <- pmin(coh, 1)  # Numerical precision

      coherence_vals[[pair_idx]] <- coh

      # Phase lag (mean phase difference, weighted by coherence)
      phase <- Arg(csd_pos)
      phase_lags[pair_idx] <- mean(phase[coh > 0.5], na.rm = TRUE)

      # Mean CSD magnitude
      csd_mags[pair_idx] <- mean(Mod(csd_pos), na.rm = TRUE)

      # Correlation of power spectra (using spectrum() output)
      if (!is.null(psds[[i]]) && !is.null(psds[[j]])) {
        spectral_corrs[pair_idx] <- tryCatch({
          cor(psds[[i]]$spec, psds[[j]]$spec, use = "complete.obs")
        }, error = function(e) NA)
      } else {
        spectral_corrs[pair_idx] <- NA
      }

      pair_idx <- pair_idx + 1
    }
  }

  # Aggregate coherence across all pairs and frequencies
  all_coherence <- unlist(coherence_vals)

  # Feature 1: Mean coherence
  coherence_mean <- mean(all_coherence, na.rm = TRUE)

  # Feature 2: Median coherence
  coherence_median <- median(all_coherence, na.rm = TRUE)

  # Feature 3: Maximum coherence
  coherence_max <- max(all_coherence, na.rm = TRUE)

  # Feature 4: Std of coherence
  coherence_std <- sd(all_coherence, na.rm = TRUE)

  # Feature 5: Fraction of high coherence (> 0.5)
  coherence_high_frac <- mean(all_coherence > 0.5, na.rm = TRUE)

  # Feature 6: Mean phase lag (absolute)
  phase_lag_mean <- mean(abs(phase_lags), na.rm = TRUE)

  # Feature 7: Std of phase lags
  phase_lag_std <- sd(phase_lags, na.rm = TRUE)

  # Feature 8: Mean spectral correlation
  spectral_corr_mean <- mean(spectral_corrs, na.rm = TRUE)

  # Feature 9: Max spectral correlation
  spectral_corr_max <- max(spectral_corrs, na.rm = TRUE)

  # Feature 10: Mean cross-spectral density magnitude
  csd_magnitude_mean <- mean(csd_mags, na.rm = TRUE)

  # Return named list
  list(
    coherence_mean = coherence_mean,
    coherence_median = coherence_median,
    coherence_max = coherence_max,
    coherence_std = coherence_std,
    coherence_high_frac = coherence_high_frac,
    phase_lag_mean = phase_lag_mean,
    phase_lag_std = phase_lag_std,
    spectral_corr_mean = spectral_corr_mean,
    spectral_corr_max = spectral_corr_max,
    csd_magnitude_mean = csd_magnitude_mean
  )
}
