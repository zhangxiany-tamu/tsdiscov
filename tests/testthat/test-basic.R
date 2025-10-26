test_that("Basic statistical features work", {
  x <- rnorm(100, mean = 5, sd = 2)

  # Test mean
  expect_equal(cpp_mean(x), mean(x), tolerance = 1e-10)

  # Test std
  expect_equal(cpp_std(x), sd(x), tolerance = 1e-10)

  # Test quantiles
  expect_equal(cpp_quantile(x, 0.5), median(x), tolerance = 1e-10)
})

test_that("ACF features work", {
  x <- rnorm(100)

  acf_vals <- cpp_acf(x, 10, TRUE)
  expect_equal(length(acf_vals), 11)
  expect_equal(acf_vals[1], 1.0, tolerance = 1e-10)
})

test_that("Entropy features work", {
  x <- rnorm(100)

  # Sample entropy should return a finite value
  se <- cpp_sample_entropy(x, 2, 0.2)
  expect_true(is.finite(se) || is.na(se))

  # Shannon entropy should be positive
  shannon <- cpp_shannon_entropy(x, 10)
  expect_true(shannon >= 0)
})

test_that("ts_features works", {
  x <- rnorm(100)

  features <- ts_features(x)
  expect_type(features, "list")
  expect_true(length(features) > 0)
  expect_true("mean" %in% names(features))
})
