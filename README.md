# tsdiscov

Time series feature extraction and selection for R with efficient C++ implementations.

## Installation

```r
# Install from GitHub
devtools::install_github("zhangxiany-tamu/tsdiscov")
```

## Quick Start

```r
library(tsdiscov)

# Extract all univariate features from a time series
x <- rnorm(200)
features <- ts_features(x, features = "all")
length(features)  # Up to 352 features (some may be NA for short/unsuitable series)

# Extract specific feature sets
basic <- ts_features(x, features = "basic")
entropy <- ts_features(x, features = "entropy")

# Bivariate features (analyze relationships between two series)
x <- rnorm(200)
y <- rnorm(200)
ccf_features <- ts_ccf_summary(x, y)        # Cross-correlation
coh_features <- ts_coherence_simple(x, y)   # Coherence analysis
granger <- ts_granger_lite(x, y)            # Granger causality
```

## Examples

### Time Series Classification

```r
# Distinguish between different types of processes
set.seed(123)

# Class 0: White noise (stationary, no structure)
white_noise <- lapply(1:30, function(i) rnorm(200))

# Class 1: Random walks (non-stationary, strong trend)
random_walks <- lapply(1:30, function(i) cumsum(rnorm(200)))

# Extract features
X <- rbind(
  do.call(rbind, lapply(white_noise, ts_features)),
  do.call(rbind, lapply(random_walks, ts_features))
)
y <- c(rep(0, 30), rep(1, 30))

# Select discriminative features
result <- suppressWarnings(select_features(X, y, fdr_level = 0.1))
cat("Selected", result$n_features_selected, "out of", result$n_features_original, "features\n")

# Use selected features for classification
library(randomForest)
model <- randomForest(result$X_selected, as.factor(y))
```

### Real-World Examples

```r
# Extract features from real-world time series
library(tsdl)

# 1. Economic series: Iowa nonfarm income
iowa_income <- tsdl[[1]]  # Quarterly, n=128
economic_features <- ts_features(iowa_income)
economic_features[c("mean", "trend_slope", "seasonal_strength", "acf_lag1")]

# 2. Agricultural: Straw yields
straw_yields <- tsdl[[50]]  # Annual, n=74
agri_features <- ts_features(straw_yields)
agri_features[c("mean", "std", "spectral_entropy", "hurst_exponent")]

# 3. Industrial production: Blooms production
blooms_prod <- tsdl[[100]]  # Quarterly, n=147
prod_features <- ts_features(blooms_prod)
prod_features[c("mean", "trend_slope", "seasonal_strength", "hurst_exponent")]

# 4. Transport: US air passenger miles
air_miles <- tsdl[[200]]  # Monthly, n=216
transport_features <- ts_features(air_miles)
transport_features[c("mean", "trend_slope", "seasonal_strength", "spectral_entropy")]

# 5. Compare multiple series from different domains
ts_list <- list(
  economic = iowa_income,
  agricultural = straw_yields,
  industrial = blooms_prod,
  transport = air_miles
)

feature_matrix <- do.call(rbind, lapply(ts_list, ts_features))
rownames(feature_matrix) <- names(ts_list)

# Compare key characteristics across domains
feature_matrix[, c("mean", "std", "trend_slope", "seasonal_strength", "spectral_entropy")]
```

## Feature Selection

Select statistically significant features using the FRESH algorithm with FDR correction:

```r
set.seed(123)

# Generate multiple time series for each class
# Class 0: White noise (20 samples)
class0 <- lapply(1:20, function(i) rnorm(150))

# Class 1: AR process (20 samples)
class1 <- lapply(1:20, function(i) as.numeric(arima.sim(list(ar = 0.8), 150)))

# Extract features and create feature matrix
X <- rbind(
  do.call(rbind, lapply(class0, ts_features)),
  do.call(rbind, lapply(class1, ts_features))
)
y <- c(rep(0, 20), rep(1, 20))  # Class labels

# Select relevant features
result <- select_features(X, y, fdr_level = 0.05)
cat(sprintf("Selected %d out of %d features\n",
            result$n_features_selected, result$n_features_original))

# Access results
selected <- result$X_selected      # Reduced feature matrix (40 x 185)
relevance <- result$relevance_table  # Statistical test results

# View top discriminative features
head(relevance[order(relevance$p_value), c("feature", "p_value", "type")], 10)
```

The feature selection implementation includes:
- Automatic selection of appropriate statistical tests (Fisher's exact, Mann-Whitney U, Kendall's tau, KS test)
- Multiple testing correction: FDR (Benjamini-Hochberg, Benjamini-Yekutieli) and FWER (Bonferroni, Holm, Hochberg, Hommel)
- Support for binary classification, multiclass classification, and regression
- Integration with standard R modeling workflows

## Features

The package provides **352 univariate features** and **6 bivariate features** across multiple categories:

### Univariate Features (352)

**Core Statistics** (20+ features): Mean, standard deviation, skewness, kurtosis, quantiles, trimmed statistics, tail asymmetry, robust shape measures

**Autocorrelation** (15+ features): ACF and PACF at various lags, timescale, first minimum, sum of squares, integral measures

**Trend Analysis** (12+ features): Linear trend parameters, mean change, second derivative, robust trend measures, monotonicity

**Frequency Domain** (40+ features): FFT coefficients, Welch power spectral density, spectral entropy, spectral shape, spectral moments, harmonic analysis, peak structure

**Scaling & Stationarity** (25+ features): Hurst exponent, stability, lumpiness, STL decomposition, differencing characteristics, unit root tests (KPSS, Phillips-Perron, ADF), rolling stability

**Complexity & Entropy** (30+ features): Sample entropy, approximate entropy, permutation entropy, Lempel-Ziv complexity, Shannon entropy, Renyi entropy, dispersion entropy, multiscale measures

**Canonical Features** (22 features): Research-validated canonical time series characteristics for robust analysis

**Change Point Detection** (8+ features): Pettitt test, CUSUM, Page-Hinkley test, segmentation analysis

**Pattern Discovery** (15+ features): Matrix profile, motif detection, recurrence analysis, Poincare maps, Friedrich coefficients

**Data Quality** (8+ features): Missing data gaps, constant segments, outlier detection, range statistics

**Nonlinear Dynamics** (20+ features): Langevin coefficients, BDS test, time reversal asymmetry, fractal dimensions (Higuchi, Katz, Petrosian), wavelet analysis

**Additional Features** (135+ features): Including returns analysis, volatility clustering, dwell times, quantile spreads, shape factors, extrema analysis, and more

### Bivariate Features (6)

**Cross-Correlation** (2 features): Maximum absolute CCF, lag at maximum

**Coherence Analysis** (2 features): Peak coherence, mean coherence across frequencies

**Granger Causality** (2 features): AIC difference, F-test p-value for predictive relationships

All computationally intensive operations are implemented in C++ via Rcpp for efficiency.

## Feature Sets

Extract features by category for targeted analysis:

```r
# All features
all_features <- ts_features(x)  # All 352 features

# Core Statistics & Distribution
ts_features(x, features = "stats")           # 7 basic statistics
ts_features(x, features = "trimmed_stats")   # 3 trimmed mean/variance
ts_features(x, features = "winsor_stats")    # 3 winsorized statistics
ts_features(x, features = "robust_shape")    # 3 robust shape measures
ts_features(x, features = "tail_ratios")     # 3 tail distribution measures
ts_features(x, features = "iqr_props")       # 4 IQR-based features
ts_features(x, features = "gini")            # 1 Gini coefficient

# Autocorrelation & Dependencies
ts_features(x, features = "acf")             # 6 ACF features
ts_features(x, features = "pacf")            # 6 PACF features
ts_features(x, features = "pacf_extended")   # 5 extended PACF analysis
ts_features(x, features = "acf_integrals")   # 2 ACF integral measures
ts_features(x, features = "acf_summary")     # 4 ACF summary statistics
ts_features(x, features = "ar")              # 3 AR model features
ts_features(x, features = "ar1_halflife")    # 1 AR(1) half-life

# Trend & Change Analysis
ts_features(x, features = "trend")           # 7 linear trend features
ts_features(x, features = "trend_extended")  # 5 extended trend measures
ts_features(x, features = "robust_trend")    # 3 robust trend estimation
ts_features(x, features = "diff")            # 6 differencing features
ts_features(x, features = "diff_stats")      # 7 differenced series stats
ts_features(x, features = "cusum")           # 2 CUSUM statistics
ts_features(x, features = "cusum_sq")        # 2 squared CUSUM
ts_features(x, features = "pettitt")         # 3 Pettitt change point test
ts_features(x, features = "page_hinkley")    # 2 Page-Hinkley test
ts_features(x, features = "monotonicity")    # 4 monotonicity measures

# Frequency Domain & Spectral Analysis
ts_features(x, features = "fft")             # 15 FFT coefficients
ts_features(x, features = "welch")           # 4 Welch PSD features
ts_features(x, features = "spectral_shape")  # 6 spectral shape measures
ts_features(x, features = "spectral_moments") # 4 spectral moments
ts_features(x, features = "spectral_psr")    # 1 power spectral ratio
ts_features(x, features = "spectral_bands")  # 5 frequency band powers
ts_features(x, features = "spectral_crest")  # 1 spectral crest factor
ts_features(x, features = "spectral_qfactor") # 1 Q-factor
ts_features(x, features = "spectral_slope_r2") # 1 spectral slope fit
ts_features(x, features = "harmonicity_index") # 1 harmonicity measure
ts_features(x, features = "fourier1")        # 4 first Fourier component
ts_features(x, features = "fourier1_resid")  # 3 Fourier residual stats

# Entropy & Complexity
ts_features(x, features = "entropy")         # 4 basic entropy measures
ts_features(x, features = "complexity")      # 10 C3 complexity features
ts_features(x, features = "renyi2_entropy")  # 1 Renyi entropy (order 2)
ts_features(x, features = "permutation_entropy") # 3 permutation entropy
ts_features(x, features = "multiscale_entropy")  # 3 multiscale entropy
ts_features(x, features = "entropy_dispersion")  # 2 dispersion entropy
ts_features(x, features = "entropy_rate_order1") # 1 entropy rate
ts_features(x, features = "multiscale_perm") # 4 multiscale permutation

# Stationarity & Scaling
ts_features(x, features = "scaling")         # 3 Hurst/stability/lumpiness
ts_features(x, features = "stl")             # 18 STL decomposition
ts_features(x, features = "stl_simple")      # 2 simple STL features
ts_features(x, features = "adf")             # 3 ADF stationarity test
ts_features(x, features = "heterogeneity")   # 3 heterogeneity measures
ts_features(x, features = "rolling_stability") # 2 rolling window stability

# Structural Features
ts_features(x, features = "structure")       # 14 structural characteristics
ts_features(x, features = "extrema")         # 5 extrema detection
ts_features(x, features = "runs")            # 4 runs test statistics
ts_features(x, features = "turning_points")  # 2 turning point counts
ts_features(x, features = "peak_intervals")  # 3 peak interval statistics
ts_features(x, features = "peak_symmetry")   # 3 peak symmetry measures
ts_features(x, features = "dwell_times")     # 4 state dwell times
ts_features(x, features = "crossings_mean_median") # 2 level crossings

# Pattern Discovery & Motifs
ts_features(x, features = "matrixprofile")   # 6 matrix profile features
ts_features(x, features = "recurrence")      # 4 recurrence analysis
ts_features(x, features = "poincare")        # 3 Poincare map features

# Nonlinear & Fractal Analysis
ts_features(x, features = "friedrich")       # 3 Friedrich coefficients
ts_features(x, features = "langevin")        # 3 Langevin equation
ts_features(x, features = "fractal_katz")    # 1 Katz fractal dimension
ts_features(x, features = "fractal_petrosian") # 1 Petrosian FD
ts_features(x, features = "fractal_higuchi") # 1 Higuchi FD
ts_features(x, features = "bds")             # 2 BDS test for independence
ts_features(x, features = "time_reversal_multi") # 3 time reversal asymmetry

# Statistical Tests
ts_features(x, features = "stattests")       # 4 statistical tests
ts_features(x, features = "ljungbox")        # 2 Ljung-Box test
ts_features(x, features = "goodness_normality") # 2 normality tests

# Miscellaneous & Specialized
ts_features(x, features = "misc")            # 3 miscellaneous features
ts_features(x, features = "shifts")          # 6 level/variance shift detection
ts_features(x, features = "holt")            # 3 Holt-Winters features
ts_features(x, features = "aggregations")    # 11 simple aggregations
ts_features(x, features = "cwt")             # 4 continuous wavelet transform
ts_features(x, features = "basic")           # 3 basic time series properties
ts_features(x, features = "frequency")       # 1 dominant frequency
ts_features(x, features = "forecastability") # 1 forecastability measure
ts_features(x, features = "arima_diag")      # 3 ARIMA diagnostic features
ts_features(x, features = "hjorth")          # 3 Hjorth parameters
ts_features(x, features = "line_length")     # 1 curve length
ts_features(x, features = "ssc")             # 1 slope sign changes
ts_features(x, features = "qcd")             # 2 quantile crossing density
ts_features(x, features = "amdf")            # 1 average magnitude difference
ts_features(x, features = "tkao")            # 2 Teager-Kaiser energy
ts_features(x, features = "tail_index")      # 2 tail index estimation
ts_features(x, features = "tail_asymmetry")  # 3 tail asymmetry measures
ts_features(x, features = "volatility_cluster") # 3 volatility clustering
ts_features(x, features = "zcr_diff")        # 1 zero-crossing rate
ts_features(x, features = "iat")             # 3 inter-arrival times
ts_features(x, features = "return_stats")    # 5 financial returns statistics
ts_features(x, features = "na_gaps")         # 3 missing data patterns
ts_features(x, features = "const_segments")  # 2 constant segment detection
ts_features(x, features = "value_range")     # 2 value range features
ts_features(x, features = "shape_factors")   # 4 shape factor measures
ts_features(x, features = "multipeak_psd")   # 3 multi-peak PSD analysis
ts_features(x, features = "block_var_ratio") # 1 block variance ratio
ts_features(x, features = "autocov_zero_lag_norm") # 1 normalized autocovariance
ts_features(x, features = "event_rate_over_threshold") # 2 event rate features
ts_features(x, features = "gph_d")           # 2 GPH estimator
ts_features(x, features = "mase")            # 1 MASE error metric
ts_features(x, features = "quantile_spread_trend") # 2 quantile spread trends

# Canonical Feature Sets (research-validated)
ts_features(x, features = "catch22")         # 20 canonical features
ts_features(x, features = "compengine")      # 12 CompEngine features

# Combine multiple feature sets
ts_features(x, features = c("stats", "acf", "entropy"))  # Custom combination
```

## Performance

Core computations use optimized C++ implementations. No external system dependencies required beyond standard R packages.

## Use Cases

- Time series classification and clustering
- Feature engineering for forecasting models
- Exploratory time series analysis
- Automated feature selection for machine learning
- Signal processing and pattern recognition
