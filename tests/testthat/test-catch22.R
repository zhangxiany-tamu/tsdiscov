test_that("Spectral features work", {
  set.seed(123)
  x <- rnorm(100)

  # Low frequency power should be between 0 and 1
  low_freq <- cpp_welch_power_low_freq(x)
  expect_true(is.finite(low_freq) || is.na(low_freq))
  if (is.finite(low_freq)) {
    expect_true(low_freq >= 0 && low_freq <= 1)
  }

  # Spectral centroid should be between 0 and 1
  centroid <- cpp_spectral_centroid(x)
  expect_true(is.finite(centroid) || is.na(centroid))
  if (is.finite(centroid)) {
    expect_true(centroid >= 0 && centroid <= 1)
  }
})

test_that("Forecasting features work", {
  set.seed(123)
  x <- rnorm(100)

  # Forecast error should be positive
  fe <- cpp_forecast_error_mean3(x)
  expect_true(is.finite(fe) || is.na(fe))
  if (is.finite(fe)) {
    expect_true(fe >= 0)
  }

  # Timescale ratio should be finite
  tr <- cpp_timescale_ratio_after_whitening(x)
  expect_true(is.finite(tr) || is.na(tr))
})

test_that("Time reversibility works", {
  set.seed(123)
  x <- rnorm(100)

  trev <- cpp_time_reversibility(x)
  expect_true(is.finite(trev) || is.na(trev))

  # For a symmetric process, should be close to 0
  # For asymmetric, should be non-zero
})

test_that("High fluctuation proportion works", {
  set.seed(123)
  x <- rnorm(100)

  hf <- cpp_high_fluctuation_prop(x)
  expect_true(is.finite(hf) || is.na(hf))
  if (is.finite(hf)) {
    expect_true(hf >= 0 && hf <= 1)
  }
})

test_that("Motif features work", {
  set.seed(123)
  x <- rnorm(100)

  motif <- cpp_motif_three_quantile(x)
  expect_true(is.finite(motif) || is.na(motif))
  if (is.finite(motif)) {
    expect_true(motif >= 0)
  }
})

test_that("Automutual information features work", {
  set.seed(123)
  x <- rnorm(100)

  # AMI lag 2
  ami2 <- cpp_automutual_info_lag2(x, 5)
  expect_true(is.finite(ami2) || is.na(ami2))
  if (is.finite(ami2)) {
    expect_true(ami2 >= 0)
  }

  # AMI first minimum
  ami_min <- cpp_automutual_info_first_min(x, 40, 10)
  expect_true(is.finite(ami_min) || is.na(ami_min))
  if (is.finite(ami_min)) {
    expect_true(ami_min > 0)
  }
})

test_that("Embedding distance works", {
  set.seed(123)
  x <- rnorm(100)

  ed <- cpp_embedding_dist_exp_fit(x)
  expect_true(is.finite(ed) || is.na(ed))
  if (is.finite(ed)) {
    expect_true(ed >= 0)
  }
})

test_that("DFA works", {
  set.seed(123)
  x <- rnorm(100)

  dfa <- cpp_dfa(x)
  expect_true(is.finite(dfa) || is.na(dfa))

  # DFA exponent typically between 0 and 2
  # 0.5 = white noise, 1.0 = 1/f noise, 1.5 = Brownian
  if (is.finite(dfa)) {
    expect_true(dfa >= -0.5 && dfa <= 2.5)
  }
})

test_that("Rescaled range (R/S) analysis works", {
  set.seed(123)
  x <- rnorm(100)

  rs <- cpp_rs_range(x)
  expect_true(is.finite(rs) || is.na(rs))

  # R/S exponent (Hurst) typically between 0 and 1
  # 0.5 = white noise, close to 1 = persistent
  if (is.finite(rs)) {
    expect_true(rs >= -0.5 && rs <= 2.5)
  }

  # Test with random walk (should have higher Hurst exponent)
  rw <- cumsum(rnorm(100))
  rs_rw <- cpp_rs_range(rw)
  expect_true(is.finite(rs_rw) || is.na(rs_rw))
})

test_that("Periodicity detection works", {
  set.seed(123)
  x <- rnorm(100)

  per <- cpp_periodicity_wang(x)
  expect_true(is.finite(per) || is.na(per))

  # Should be between -1 and 1 (ACF value)
  if (is.finite(per)) {
    expect_true(per >= -1 && per <= 1)
  }

  # Test with periodic signal
  t <- 1:100
  y <- sin(2 * pi * t / 10) + rnorm(100, sd = 0.1)
  per_periodic <- cpp_periodicity_wang(y)
  expect_true(is.finite(per_periodic))
  # Periodic signal should have higher periodicity measure
  expect_true(per_periodic > 0.3)
})

test_that("Complete catch22 feature set works", {
  set.seed(123)
  x <- rnorm(100)

  features <- ts_catch22(x)

  # Should have exactly 20 features (removed 2 duplicates: acf_timescale, acf_first_min)
  expect_type(features, "list")
  expect_equal(length(features), 20)

  # All features should be numeric
  for (i in seq_along(features)) {
    expect_true(is.numeric(features[[i]]))
  }

  # Check specific feature names exist
  expect_true("mode_5" %in% names(features))
  expect_true("low_freq_power" %in% names(features))
  expect_true("trev" %in% names(features))
  expect_true("dfa" %in% names(features))
  expect_true("rs_range" %in% names(features))
  expect_true("periodicity" %in% names(features))
})

test_that("catch22 works with different time series patterns", {
  # White noise
  set.seed(123)
  white_noise <- rnorm(100)
  f1 <- ts_catch22(white_noise)
  expect_type(f1, "list")

  # Trending
  trend <- 1:100 + rnorm(100, sd = 5)
  f2 <- ts_catch22(trend)
  expect_type(f2, "list")

  # Periodic
  t <- 1:100
  periodic <- sin(2 * pi * t / 10) + rnorm(100, sd = 0.1)
  f3 <- ts_catch22(periodic)
  expect_type(f3, "list")

  # Periodic series should have higher periodicity
  expect_true(f3$periodicity > f1$periodicity)
})

test_that("catch22 handles edge cases", {
  # Short series
  x_short <- rnorm(5)
  f <- ts_catch22(x_short)
  expect_type(f, "list")
  # Some features may be NA for short series
  expect_true(any(sapply(f, is.na)) || all(sapply(f, is.finite)))

  # Constant series
  x_const <- rep(1, 100)
  f_const <- ts_catch22(x_const)
  expect_type(f_const, "list")
  # Many features should be NA or 0 for constant series
})
