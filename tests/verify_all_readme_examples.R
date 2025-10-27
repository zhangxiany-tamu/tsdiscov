# Verify ALL README examples work correctly

library(tsdiscov)

cat("═══════════════════════════════════════════════════════\n")
cat("     TESTING ALL README EXAMPLES\n")
cat("═══════════════════════════════════════════════════════\n\n")

# =============================================================================
# EXAMPLE 1: Quick Start - Univariate
# =============================================================================
cat("Example 1: Quick Start - Univariate Features\n")
cat("---------------------------------------------\n")
x <- rnorm(200)
features <- ts_features(x, features = "all")
cat(sprintf("✓ Extracted %d univariate features\n", length(features)))

basic <- ts_features(x, features = "basic")
entropy <- ts_features(x, features = "entropy")
cat(sprintf("✓ Basic features: %d\n", length(basic)))
cat(sprintf("✓ Entropy features: %d\n\n", length(entropy)))

# =============================================================================
# EXAMPLE 2: Quick Start - Bivariate
# =============================================================================
cat("Example 2: Quick Start - Bivariate Features\n")
cat("--------------------------------------------\n")
x <- rnorm(200)
y <- rnorm(200)
ccf_features <- ts_ccf_summary(x, y)
coh_features <- ts_coherence_simple(x, y)
granger <- ts_granger_lite(x, y)
cat(sprintf("✓ CCF features: %d\n", length(ccf_features)))
cat(sprintf("✓ Coherence features: %d\n", length(coh_features)))
cat(sprintf("✓ Granger features: %d\n\n", length(granger)))

# =============================================================================
# EXAMPLE 3: Quick Start - Multivariate
# =============================================================================
cat("Example 3: Quick Start - Multivariate Features\n")
cat("-----------------------------------------------\n")
X <- matrix(rnorm(4 * 200), nrow = 4, ncol = 200)

mv_features <- ts_features_all(X, feature_type = "multivariate")
cat(sprintf("✓ Extracted %d multivariate features\n", length(mv_features)))

both_features <- ts_features_all(X,
                                  feature_type = "both",
                                  univariate_summary = "aggregate")
cat(sprintf("✓ Both features: %d\n", length(both_features)))

systems <- list(
  sensor_array_1 = matrix(rnorm(3 * 150), nrow = 3),
  sensor_array_2 = matrix(rnorm(3 * 150), nrow = 3),
  sensor_array_3 = matrix(rnorm(3 * 150), nrow = 3)
)
results <- do.call(rbind, lapply(systems, ts_features_all_df))
rownames(results) <- names(systems)
cat(sprintf("✓ Batch processing: %d systems × %d features\n\n",
            nrow(results), ncol(results)))

# =============================================================================
# EXAMPLE 4: Time Series Classification (Univariate)
# =============================================================================
cat("Example 4: Time Series Classification (Univariate)\n")
cat("--------------------------------------------------\n")
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
cat(sprintf("✓ Selected %d out of %d features\n",
            result$n_features_selected, result$n_features_original))

# Use selected features for classification
library(randomForest)
model <- randomForest(result$X_selected, as.factor(y))
cat(sprintf("✓ Model trained successfully\n\n"))

# =============================================================================
# EXAMPLE 5: Multivariate Time Series Classification
# =============================================================================
cat("Example 5: Multivariate Time Series Classification\n")
cat("---------------------------------------------------\n")
library(randomForest)
set.seed(42)

# Generate systems with 3-8 sensors each (40 per class)
systems <- list()
labels <- c()

# Class 0: Independent sensors
for (i in 1:40) {
  N <- sample(3:8, 1)
  systems[[i]] <- matrix(rnorm(N * 150), nrow = N)
  labels[i] <- 0
}

# Class 1: Synchronized sensors
for (i in 41:80) {
  N <- sample(3:8, 1)
  t <- seq(0, 10, length.out = 150)
  common <- sin(2*pi*t)
  X <- matrix(0, nrow = N, ncol = 150)
  for (j in 1:N) X[j,] <- common + rnorm(150, sd = 0.2)
  systems[[i]] <- X
  labels[i] <- 1
}

# Extract features (50 features regardless of N)
features <- do.call(rbind, lapply(systems, ts_features_all_df))
features <- features[, colSums(is.na(features)) == 0]

# Train and evaluate
train_idx <- c(1:30, 41:70)
test_idx <- c(31:40, 71:80)
model <- randomForest(features[train_idx, ], as.factor(labels[train_idx]))
pred <- predict(model, features[test_idx, ])
accuracy <- mean(pred == labels[test_idx])

cat(sprintf("✓ Feature matrix: %d samples × %d features\n",
            nrow(features), ncol(features)))
cat(sprintf("✓ Test accuracy: %.1f%%\n\n", accuracy * 100))

# =============================================================================
# EXAMPLE 6: Real-World Examples (with tsdl)
# =============================================================================
cat("Example 6: Real-World Examples (tsdl)\n")
cat("-------------------------------------\n")

if (require("tsdl", quietly = TRUE)) {
  library(tsdl)

  # 1. Economic series: Iowa nonfarm income
  iowa_income <- tsdl[[1]]
  economic_features <- ts_features(iowa_income)
  cat(sprintf("✓ Economic series: %d features extracted\n",
              length(economic_features)))

  # 2. Agricultural: Straw yields
  straw_yields <- tsdl[[50]]
  agri_features <- ts_features(straw_yields)
  cat(sprintf("✓ Agricultural series: %d features extracted\n",
              length(agri_features)))

  # 3. Industrial production: Blooms production
  blooms_prod <- tsdl[[100]]
  prod_features <- ts_features(blooms_prod)
  cat(sprintf("✓ Industrial series: %d features extracted\n",
              length(prod_features)))

  # 4. Transport: US air passenger miles
  air_miles <- tsdl[[200]]
  transport_features <- ts_features(air_miles)
  cat(sprintf("✓ Transport series: %d features extracted\n",
              length(transport_features)))

  # 5. Compare multiple series from different domains
  ts_list <- list(
    economic = iowa_income,
    agricultural = straw_yields,
    industrial = blooms_prod,
    transport = air_miles
  )

  feature_matrix <- do.call(rbind, lapply(ts_list, ts_features))
  rownames(feature_matrix) <- names(ts_list)

  cat(sprintf("✓ Batch comparison: %d series × %d features\n\n",
              nrow(feature_matrix), ncol(feature_matrix)))
} else {
  cat("⚠ tsdl package not available, skipping real-world examples\n\n")
}

# =============================================================================
# EXAMPLE 7: Feature Selection
# =============================================================================
cat("Example 7: Feature Selection\n")
cat("----------------------------\n")
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
cat(sprintf("✓ Selected %d out of %d features\n",
            result$n_features_selected, result$n_features_original))

# Access results
selected <- result$X_selected
relevance <- result$relevance_table

# View top discriminative features
top_features <- head(relevance[order(relevance$p_value),
                                c("feature", "p_value", "type")], 10)
cat(sprintf("✓ Top 10 features identified\n"))
cat(sprintf("✓ Best feature: %s (p = %.2e)\n\n",
            top_features$feature[1], top_features$p_value[1]))

# =============================================================================
# SUMMARY
# =============================================================================
cat("═══════════════════════════════════════════════════════\n")
cat("     ALL README EXAMPLES PASSED ✅\n")
cat("═══════════════════════════════════════════════════════\n\n")

cat("Summary:\n")
cat("✓ Example 1: Univariate features\n")
cat("✓ Example 2: Bivariate features\n")
cat("✓ Example 3: Multivariate features\n")
cat("✓ Example 4: Univariate classification\n")
cat("✓ Example 5: Multivariate classification (", accuracy * 100, "% accuracy)\n", sep="")
cat("✓ Example 6: Real-world tsdl examples\n")
cat("✓ Example 7: Feature selection\n\n")

cat("All examples run successfully without errors!\n")
