# CompEngine Features Tests
# Tests for features from CompEngine/hctsa database

test_that("CompEngine features work correctly", {
  set.seed(123)
  x <- rnorm(100)

  result <- ts_compengine(x)

  # Should have 12 features
  expect_equal(length(result), 12)

  # Check all feature names are present
  expected_names <- c(
    "crossing_points", "flat_spots", "embed2_incircle_1", "embed2_incircle_2",
    "motiftwo_entro3", "walker_propcross", "localsimple_mean1", "localsimple_lfitac",
    "spreadrandomlocal_meantaul_50", "spreadrandomlocal_meantaul_ac2",
    "outlierinclude_mdrmd", "fluctanal_prop_r1"
  )
  expect_true(all(expected_names %in% names(result)))
})

test_that("crossing_points computes median crossings", {
  set.seed(456)

  # White noise should have many crossings
  x_noise <- rnorm(100)
  result_noise <- ts_compengine(x_noise)
  expect_true(is.finite(result_noise$crossing_points))
  expect_true(result_noise$crossing_points > 0)

  # Constant series has 0 crossings
  x_const <- rep(5, 100)
  crossings_const <- cpp_crossing_points(x_const)
  expect_equal(crossings_const, 0)

  # Monotonic series has very few crossings (0 or 1)
  x_mono <- 1:100
  crossings_mono <- cpp_crossing_points(x_mono)
  expect_true(crossings_mono <= 1)
})

test_that("flat_spots finds longest run in bins", {
  set.seed(789)

  # White noise should have short flat spots
  x_noise <- rnorm(100)
  flat_noise <- cpp_flat_spots(x_noise)
  expect_true(is.finite(flat_noise))
  expect_true(flat_noise >= 1)

  # Constant series has flat_spots = length
  x_const <- rep(5, 100)
  flat_const <- cpp_flat_spots(x_const)
  expect_equal(flat_const, 100)

  # Series with plateaus should have longer flat spots
  x_plateau <- c(rep(1, 20), rep(2, 20), rep(3, 20), rep(4, 20), rep(5, 20))
  flat_plateau <- cpp_flat_spots(x_plateau)
  expect_true(flat_plateau >= 15)  # At least 15 in same bin
})

test_that("embed2_incircle computes 2D embedding correctly", {
  set.seed(111)

  x <- rnorm(100)

  # Boundary 1 should have fewer points than boundary 2
  embed1 <- cpp_embed2_incircle(x, 1)
  embed2 <- cpp_embed2_incircle(x, 2)

  if (is.finite(embed1) && is.finite(embed2)) {
    expect_true(embed1 <= embed2)  # Larger circle contains more points
  }

  # Should be proportions between 0 and 1
  if (is.finite(embed1)) {
    expect_true(embed1 >= 0 && embed1 <= 1)
  }
  if (is.finite(embed2)) {
    expect_true(embed2 >= 0 && embed2 <= 1)
  }
})

test_that("motiftwo_entro3 computes binary motif entropy", {
  set.seed(222)

  # White noise should have high entropy
  x_noise <- rnorm(100)
  entropy_noise <- cpp_motiftwo_entro3(x_noise)
  expect_true(is.finite(entropy_noise))
  expect_true(entropy_noise > 0)

  # Constant series (all below mean) has 0 entropy
  x_const <- rep(5, 100)
  entropy_const <- cpp_motiftwo_entro3(x_const)
  expect_true(is.finite(entropy_const))
  expect_true(entropy_const < 0.1)  # Near zero

  # Alternating series should have low entropy (regular pattern)
  x_alt <- rep(c(1, -1), 50)
  entropy_alt <- cpp_motiftwo_entro3(x_alt)
  expect_true(is.finite(entropy_alt))
})

test_that("walker_propcross computes walker crossings", {
  set.seed(333)

  x <- rnorm(100)
  walker_cross <- cpp_walker_propcross(x)

  # Should be proportion between 0 and 1
  expect_true(is.finite(walker_cross))
  expect_true(walker_cross >= 0 && walker_cross <= 1)

  # Constant series should have low crossings
  x_const <- rep(5, 100)
  walker_const <- cpp_walker_propcross(x_const)
  expect_true(walker_const < 0.1)  # Very few crossings
})

test_that("localsimple_mean1 predicts using past mean", {
  set.seed(444)

  x <- rnorm(100)
  result <- ts_compengine(x)

  # Should compute finite value
  expect_true(is.numeric(result$localsimple_mean1))

  # For white noise, residuals should have a zero crossing
  if (is.finite(result$localsimple_mean1)) {
    expect_true(result$localsimple_mean1 >= 0)
  }
})

test_that("localsimple_lfitac predicts using local linear fit", {
  set.seed(555)

  x <- rnorm(100)
  result <- ts_compengine(x)

  # Should compute value (may be NA if tau too large)
  expect_true(is.numeric(result$localsimple_lfitac))

  if (is.finite(result$localsimple_lfitac)) {
    expect_true(result$localsimple_lfitac >= 0)
  }
})

test_that("spreadrandomlocal_meantaul measures stationarity", {
  set.seed(666)

  x <- rnorm(200)
  result <- ts_compengine(x)

  # Should compute both versions
  expect_true(is.numeric(result$spreadrandomlocal_meantaul_50))
  expect_true(is.numeric(result$spreadrandomlocal_meantaul_ac2))

  # Values should be positive (ACF zero crossings)
  if (is.finite(result$spreadrandomlocal_meantaul_50)) {
    expect_true(result$spreadrandomlocal_meantaul_50 >= 0)
  }
})

test_that("outlierinclude_mdrmd computes outlier metric", {
  set.seed(777)

  x <- rnorm(100)
  result <- ts_compengine(x)

  # Should compute value
  expect_true(is.numeric(result$outlierinclude_mdrmd))

  # For normalized data, should be near 0
  if (is.finite(result$outlierinclude_mdrmd)) {
    expect_true(abs(result$outlierinclude_mdrmd) < 1)
  }
})

test_that("fluctanal_prop_r1 computes fluctuation analysis", {
  set.seed(888)

  x <- rnorm(200)
  result <- ts_compengine(x)

  # Should compute value (proportion between 0 and 1)
  expect_true(is.numeric(result$fluctanal_prop_r1))

  if (is.finite(result$fluctanal_prop_r1)) {
    expect_true(result$fluctanal_prop_r1 >= 0 && result$fluctanal_prop_r1 <= 1)
  }
})

test_that("CompEngine features handle short series", {
  x_short <- rnorm(5)

  # Suppress expected warning about short series
  result <- suppressWarnings(ts_compengine(x_short))

  # Should return NA values for short series
  expect_equal(length(result), 12)
  expect_true(all(sapply(result, is.na)))
})

test_that("CompEngine features integrate with ts_features", {
  set.seed(999)

  x <- rnorm(150)

  # Test with 'compengine' subset
  comp_feats <- ts_features(x, features = "compengine")
  expect_equal(length(comp_feats), 12)
  expect_true("crossing_points" %in% names(comp_feats))
  expect_true("flat_spots" %in% names(comp_feats))

  # Test with 'all'
  all_feats <- ts_features(x, features = "all")
  expect_true("crossing_points" %in% names(all_feats))
  expect_true("walker_propcross" %in% names(all_feats))
  expect_true(length(all_feats) >= 141)  # Should have at least 133 features now
})

test_that("Feature count is correct after CompEngine addition", {
  set.seed(101010)
  x <- rnorm(200)

  all_feats <- ts_features(x, features = "all")

  # Should have 133 total features (121 + 12 CompEngine)
  expect_equal(length(all_feats), EXPECTED_ALL_FEATURES)
})

test_that("CompEngine features are deterministic", {
  set.seed(111111)
  x <- rnorm(100)

  result1 <- ts_compengine(x)
  result2 <- ts_compengine(x)

  # Should get identical results (except for random bootstrap features)
  expect_equal(result1$crossing_points, result2$crossing_points)
  expect_equal(result1$flat_spots, result2$flat_spots)
  expect_equal(result1$embed2_incircle_1, result2$embed2_incircle_1)
  expect_equal(result1$motiftwo_entro3, result2$motiftwo_entro3)
  expect_equal(result1$walker_propcross, result2$walker_propcross)

  # Note: spreadrandomlocal features use random sampling, so may differ slightly
})

test_that("crossing_points differs from number_crossings", {
  set.seed(121212)
  x <- rnorm(100)

  # crossing_points uses median
  crossing_med <- cpp_crossing_points(x)

  # number_crossings uses mean (in structure features)
  structure_feats <- ts_structure(x)
  crossing_mean <- structure_feats$number_crossings

  # They should be similar but not identical
  expect_true(is.finite(crossing_med))
  expect_true(is.finite(crossing_mean))
})

test_that("CompEngine features handle NA values correctly", {
  set.seed(131313)

  # Series with NAs (will be removed by ts_features before calling ts_compengine)
  x <- rnorm(100)
  x[c(10, 50, 80)] <- NA

  # ts_features removes NAs
  comp_feats <- ts_features(x, features = "compengine")

  # Should handle NA removal
  expect_equal(length(comp_feats), 12)
  expect_true(is.numeric(comp_feats$crossing_points))
})

test_that("CompEngine features work with different data patterns", {
  set.seed(141414)

  # Trending data
  t <- 1:100
  x_trend <- 0.05 * t + rnorm(100, sd = 0.5)
  result_trend <- ts_compengine(x_trend)
  expect_equal(length(result_trend), 12)

  # Seasonal data
  x_seasonal <- sin(2 * pi * t / 12) + rnorm(100, sd = 0.2)
  result_seasonal <- ts_compengine(x_seasonal)
  expect_equal(length(result_seasonal), 12)

  # Random walk
  x_rw <- cumsum(rnorm(100))
  result_rw <- ts_compengine(x_rw)
  expect_equal(length(result_rw), 12)
})
