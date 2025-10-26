# Test helper: Expected feature counts
#
# This file defines constants for expected feature counts.
# When features are added to the package, update these constants by running:
#   devtools::load_all()
#   length(list_ts_feature_names(rnorm(100), "all"))
#   length(list_ts_feature_names(rnorm(100), "stats"))
#   etc.

# Expected count for "all" features (352 as of 2025-01-26)
EXPECTED_ALL_FEATURES <- 352

# Expected counts for specific feature sets (commonly tested)
EXPECTED_STATS_FEATURES <- 7
EXPECTED_ACF_FEATURES <- 6
EXPECTED_MISC_FEATURES <- 3
EXPECTED_SHIFTS_FEATURES <- 6
