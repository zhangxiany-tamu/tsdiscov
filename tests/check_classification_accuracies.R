# Check classification accuracies for Examples 4 and 5

library(tsdiscov)
library(randomForest)

cat("═══════════════════════════════════════════════════════\n")
cat("  CLASSIFICATION ACCURACIES - README EXAMPLES\n")
cat("═══════════════════════════════════════════════════════\n\n")

# =============================================================================
# EXAMPLE 4: Univariate Time Series Classification
# =============================================================================
cat("Example 4: Univariate Time Series Classification\n")
cat("-------------------------------------------------\n")
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
cat(sprintf("Selected features: %d out of %d\n",
            result$n_features_selected, result$n_features_original))

# Train model with selected features
model <- randomForest(result$X_selected, as.factor(y))

# Get predictions
pred <- predict(model, result$X_selected)
accuracy <- mean(pred == y)

cat(sprintf("\nTraining Accuracy: %.1f%%\n", accuracy * 100))
cat(sprintf("Confusion Matrix:\n"))
cm <- table(Predicted = pred, Actual = y)
print(cm)
cat("\n")

# =============================================================================
# EXAMPLE 5: Multivariate Time Series Classification
# =============================================================================
cat("Example 5: Multivariate Time Series Classification\n")
cat("---------------------------------------------------\n")
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

# Extract features
features <- do.call(rbind, lapply(systems, ts_features_all_df))
features <- features[, colSums(is.na(features)) == 0]

cat(sprintf("Feature matrix: %d samples × %d features\n",
            nrow(features), ncol(features)))

# Train-test split
train_idx <- c(1:30, 41:70)
test_idx <- c(31:40, 71:80)

# Train model
model <- randomForest(features[train_idx, ], as.factor(labels[train_idx]))

# Training accuracy
pred_train <- predict(model, features[train_idx, ])
train_acc <- mean(pred_train == labels[train_idx])

# Test accuracy
pred_test <- predict(model, features[test_idx, ])
test_acc <- mean(pred_test == labels[test_idx])

cat(sprintf("\nTraining Accuracy: %.1f%%\n", train_acc * 100))
cat(sprintf("Test Accuracy: %.1f%%\n", test_acc * 100))
cat(sprintf("\nTest Confusion Matrix:\n"))
cm_test <- table(Predicted = pred_test, Actual = labels[test_idx])
print(cm_test)
cat("\n")

# =============================================================================
# SUMMARY
# =============================================================================
cat("═══════════════════════════════════════════════════════\n")
cat("  SUMMARY\n")
cat("═══════════════════════════════════════════════════════\n\n")

cat(sprintf("Example 4 (Univariate Classification):\n"))
cat(sprintf("  Training Accuracy: %.1f%%\n", accuracy * 100))
cat(sprintf("  Task: Distinguish white noise vs random walks\n\n"))

cat(sprintf("Example 5 (Multivariate Classification):\n"))
cat(sprintf("  Training Accuracy: %.1f%%\n", train_acc * 100))
cat(sprintf("  Test Accuracy: %.1f%%\n", test_acc * 100))
cat(sprintf("  Task: Distinguish independent vs synchronized sensors\n"))
cat(sprintf("  Special: Handles varying sensor counts (3-8)\n\n"))
