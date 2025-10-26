#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <set>
#include <unordered_set>
using namespace Rcpp;

// [[Rcpp::export]]
double cpp_c3(NumericVector x, int lag = 1) {
    /**
     * C3 statistic to measure non-linearity in time series
     *
     * Calculates: mean(x[i+2*lag] * x[i+lag] * x[i])
     *
     * Reference: Schreiber & Schmitz (1997)
     * "Discrimination power of measures for nonlinearity in a time series"
     * PHYSICAL REVIEW E, VOL 55, NO 5
     *
     * @param x Numeric vector (time series)
     * @param lag The lag to use (default 1)
     * @return C3 statistic value
     */
    int n = x.size();

    if (n < 2 * lag + 1) {
        return NA_REAL;
    }

    double sum = 0.0;
    int count = 0;

    for (int i = 0; i < n - 2 * lag; i++) {
        if (!NumericVector::is_na(x[i]) &&
            !NumericVector::is_na(x[i + lag]) &&
            !NumericVector::is_na(x[i + 2 * lag])) {
            sum += x[i + 2 * lag] * x[i + lag] * x[i];
            count++;
        }
    }

    if (count == 0) {
        return NA_REAL;
    }

    return sum / count;
}

// [[Rcpp::export]]
double cpp_cid_ce(NumericVector x, bool normalize = true) {
    /**
     * Complexity-Invariant Distance (CID)
     *
     * Estimates time series complexity based on consecutive differences.
     * More complex series have more peaks and valleys.
     *
     * Calculates: sqrt(sum((x[i] - x[i-1])^2))
     *
     * Reference: Batista et al. (2014)
     * "CID: an efficient complexity-invariant distance for time series"
     * Data Mining and Knowledge Discovery 28.3: 634-669
     *
     * @param x Numeric vector (time series)
     * @param normalize Should the time series be z-transformed first?
     * @return CID complexity estimate
     */
    NumericVector x_copy = clone(x);
    int n = x_copy.size();

    if (n < 2) {
        return NA_REAL;
    }

    // Normalize if requested
    if (normalize) {
        double mean = 0.0;
        for (int i = 0; i < n; i++) {
            mean += x_copy[i];
        }
        mean /= n;

        double var = 0.0;
        for (int i = 0; i < n; i++) {
            double diff = x_copy[i] - mean;
            var += diff * diff;
        }
        var /= n;
        double sd = std::sqrt(var);

        if (sd > 0) {
            for (int i = 0; i < n; i++) {
                x_copy[i] = (x_copy[i] - mean) / sd;
            }
        }
    }

    // Calculate sum of squared differences
    double sum = 0.0;
    for (int i = 1; i < n; i++) {
        if (!NumericVector::is_na(x_copy[i]) && !NumericVector::is_na(x_copy[i-1])) {
            double diff = x_copy[i] - x_copy[i-1];
            sum += diff * diff;
        }
    }

    return std::sqrt(sum);
}

// [[Rcpp::export]]
double cpp_lempel_ziv_complexity(NumericVector x, int bins = 10) {
    /**
     * Lempel-Ziv Complexity
     *
     * Complexity estimate based on Lempel-Ziv compression algorithm.
     * Counts number of unique sub-words needed to encode the binned series.
     *
     * Reference: https://github.com/Naereen/Lempel-Ziv_Complexity
     *
     * @param x Numeric vector (time series)
     * @param bins Number of bins for discretization
     * @return Normalized LZ complexity (0 to 1)
     */
    int n = x.size();

    if (n < 2 || bins < 2) {
        return NA_REAL;
    }

    // Find min and max for binning
    double min_val = x[0];
    double max_val = x[0];
    for (int i = 1; i < n; i++) {
        if (x[i] < min_val) min_val = x[i];
        if (x[i] > max_val) max_val = x[i];
    }

    if (max_val == min_val) {
        return 0.0;  // Constant series has zero complexity
    }

    // Bin the data
    std::vector<int> sequence(n);
    double bin_width = (max_val - min_val) / bins;

    for (int i = 0; i < n; i++) {
        int bin = static_cast<int>((x[i] - min_val) / bin_width);
        if (bin >= bins) bin = bins - 1;
        sequence[i] = bin;
    }

    // Lempel-Ziv algorithm
    std::set<std::vector<int>> sub_strings;
    int ind = 0;
    int inc = 1;

    while (ind + inc <= n) {
        // Get substring
        std::vector<int> sub_str(sequence.begin() + ind, sequence.begin() + ind + inc);

        if (sub_strings.find(sub_str) != sub_strings.end()) {
            // Already seen, increase length
            inc++;
        } else {
            // New substring
            sub_strings.insert(sub_str);
            ind += inc;
            inc = 1;
        }
    }

    return static_cast<double>(sub_strings.size()) / n;
}

// [[Rcpp::export]]
double cpp_index_mass_quantile(NumericVector x, double q = 0.5) {
    /**
     * Index Mass Quantile
     *
     * Calculates the relative index where q% of the mass of x lies to the left.
     * For q=0.5, this returns the mass center of the time series.
     *
     * @param x Numeric vector (time series)
     * @param q Quantile (0 to 1)
     * @return Normalized index (0 to 1) where q% of mass is to the left
     */
    int n = x.size();

    if (n == 0 || q < 0 || q > 1) {
        return NA_REAL;
    }

    // Calculate absolute values
    std::vector<double> abs_x(n);
    double total_mass = 0.0;

    for (int i = 0; i < n; i++) {
        abs_x[i] = std::abs(x[i]);
        total_mass += abs_x[i];
    }

    if (total_mass == 0) {
        return NA_REAL;
    }

    // Calculate cumulative mass
    double cumsum = 0.0;
    for (int i = 0; i < n; i++) {
        cumsum += abs_x[i];
        if (cumsum / total_mass >= q) {
            return static_cast<double>(i + 1) / n;
        }
    }

    return 1.0;
}

// [[Rcpp::export]]
double cpp_mean_n_absolute_max(NumericVector x, int number_of_maxima = 1) {
    /**
     * Mean of N Absolute Maximum Values
     *
     * Calculates the arithmetic mean of the n absolute maximum values.
     * Robust alternative to maximum that considers top-n values.
     *
     * @param x Numeric vector (time series)
     * @param number_of_maxima Number of maxima to average
     * @return Mean of n largest absolute values
     */
    int n = x.size();

    if (number_of_maxima <= 0 || n < number_of_maxima) {
        return NA_REAL;
    }

    // Create vector of absolute values
    std::vector<double> abs_values(n);
    for (int i = 0; i < n; i++) {
        abs_values[i] = std::abs(x[i]);
    }

    // Partial sort to get top number_of_maxima
    std::partial_sort(abs_values.begin(),
                      abs_values.begin() + number_of_maxima,
                      abs_values.end(),
                      std::greater<double>());

    // Calculate mean of top values
    double sum = 0.0;
    for (int i = 0; i < number_of_maxima; i++) {
        sum += abs_values[i];
    }

    return sum / number_of_maxima;
}

// [[Rcpp::export]]
double cpp_energy_ratio_by_chunks(NumericVector x, int num_segments = 10, int segment_focus = 0) {
    /**
     * Energy Ratio by Chunks
     *
     * Calculates the sum of squares of chunk i out of N chunks, expressed as
     * a ratio with the sum of squares over the whole series.
     *
     * @param x Numeric vector (time series)
     * @param num_segments Number of segments to divide series into
     * @param segment_focus Which segment to focus on (0-indexed)
     * @return Ratio of energy in focused segment to total energy
     */
    int n = x.size();

    if (num_segments <= 0 || segment_focus < 0 || segment_focus >= num_segments) {
        return NA_REAL;
    }

    // Calculate total energy
    double total_energy = 0.0;
    for (int i = 0; i < n; i++) {
        total_energy += x[i] * x[i];
    }

    if (total_energy == 0) {
        return NA_REAL;
    }

    // Determine chunk sizes (distribute remainder across first chunks)
    int base_size = n / num_segments;
    int remainder = n % num_segments;

    // Calculate start and end indices for focused segment
    int start_idx = 0;
    for (int i = 0; i < segment_focus; i++) {
        start_idx += base_size + (i < remainder ? 1 : 0);
    }

    int segment_size = base_size + (segment_focus < remainder ? 1 : 0);
    int end_idx = start_idx + segment_size;

    // Calculate energy in focused segment
    double segment_energy = 0.0;
    for (int i = start_idx; i < end_idx && i < n; i++) {
        segment_energy += x[i] * x[i];
    }

    return segment_energy / total_energy;
}

// [[Rcpp::export]]
double cpp_change_quantiles(NumericVector x, double ql = 0.0, double qh = 1.0,
                              bool isabs = false, std::string f_agg = "mean") {
    /**
     * Change Quantiles
     *
     * Fixes a corridor given by quantiles ql and qh of distribution of x.
     * Calculates the aggregated value of consecutive changes inside this corridor.
     *
     * @param x Numeric vector (time series)
     * @param ql Lower quantile of corridor (0 to 1)
     * @param qh Higher quantile of corridor (0 to 1)
     * @param isabs Should absolute differences be taken?
     * @param f_agg Aggregator function: "mean", "std", "var", "median"
     * @return Aggregated change value inside corridor
     */
    int n = x.size();

    if (n < 2 || ql >= qh) {
        return 0.0;
    }

    // Calculate differences
    std::vector<double> diffs(n - 1);
    for (int i = 0; i < n - 1; i++) {
        diffs[i] = x[i + 1] - x[i];
        if (isabs) {
            diffs[i] = std::abs(diffs[i]);
        }
    }

    // Calculate quantiles of x
    NumericVector x_sorted = clone(x);
    std::sort(x_sorted.begin(), x_sorted.end());

    int idx_low = static_cast<int>(ql * (n - 1));
    int idx_high = static_cast<int>(qh * (n - 1));

    double q_low = x_sorted[idx_low];
    double q_high = x_sorted[idx_high];

    if (q_low >= q_high) {
        return 0.0;
    }

    // Find changes that start and end inside corridor
    std::vector<double> changes_in_corridor;
    for (int i = 0; i < n - 1; i++) {
        bool start_in = (x[i] >= q_low) && (x[i] <= q_high);
        bool end_in = (x[i + 1] >= q_low) && (x[i + 1] <= q_high);

        if (start_in && end_in) {
            changes_in_corridor.push_back(diffs[i]);
        }
    }

    if (changes_in_corridor.empty()) {
        return 0.0;
    }

    // Apply aggregator
    double result = 0.0;
    int count = changes_in_corridor.size();

    if (f_agg == "mean") {
        for (double val : changes_in_corridor) {
            result += val;
        }
        result /= count;
    } else if (f_agg == "median") {
        std::sort(changes_in_corridor.begin(), changes_in_corridor.end());
        if (count % 2 == 0) {
            result = (changes_in_corridor[count/2 - 1] + changes_in_corridor[count/2]) / 2.0;
        } else {
            result = changes_in_corridor[count/2];
        }
    } else if (f_agg == "var") {
        double mean = 0.0;
        for (double val : changes_in_corridor) {
            mean += val;
        }
        mean /= count;

        for (double val : changes_in_corridor) {
            double diff = val - mean;
            result += diff * diff;
        }
        result /= count;
    } else if (f_agg == "std") {
        double mean = 0.0;
        for (double val : changes_in_corridor) {
            mean += val;
        }
        mean /= count;

        for (double val : changes_in_corridor) {
            double diff = val - mean;
            result += diff * diff;
        }
        result = std::sqrt(result / count);
    } else {
        // Default to mean
        for (double val : changes_in_corridor) {
            result += val;
        }
        result /= count;
    }

    return result;
}

// [[Rcpp::export]]
double cpp_fourier_entropy(NumericVector x, int bins = 10) {
    /**
     * Fourier Entropy
     *
     * Calculates binned entropy of the power spectral density using Welch method.
     * Measures frequency domain complexity.
     *
     * Reference: https://hackaday.io/project/707-complexity-of-a-time-series
     *
     * @param x Numeric vector (time series)
     * @param bins Number of bins for entropy calculation
     * @return Fourier entropy value
     */
    int n = x.size();

    if (n < 4 || bins < 2) {
        return NA_REAL;
    }

    // For simplicity, use full FFT (in real application, Welch method segments the data)
    // Calculate FFT power spectrum
    int fft_size = n / 2;
    std::vector<double> power(fft_size);

    // Compute power spectrum (simplified - ideally use Welch method)
    // Here we'll approximate by computing periodogram
    double mean = 0.0;
    for (int i = 0; i < n; i++) {
        mean += x[i];
    }
    mean /= n;

    // Center the data
    std::vector<double> centered(n);
    for (int i = 0; i < n; i++) {
        centered[i] = x[i] - mean;
    }

    // Compute power spectrum magnitude (simplified periodogram)
    for (int k = 0; k < fft_size; k++) {
        double real_part = 0.0;
        double imag_part = 0.0;

        for (int i = 0; i < n; i++) {
            double angle = 2.0 * M_PI * k * i / n;
            real_part += centered[i] * std::cos(angle);
            imag_part += centered[i] * std::sin(angle);
        }

        power[k] = real_part * real_part + imag_part * imag_part;
    }

    // Normalize power spectrum
    double max_power = *std::max_element(power.begin(), power.end());
    if (max_power == 0) {
        return NA_REAL;
    }

    for (int i = 0; i < fft_size; i++) {
        power[i] /= max_power;
    }

    // Bin the power spectrum and calculate entropy
    std::vector<double> bin_counts(bins, 0.0);
    double bin_width = 1.0 / bins;

    for (int i = 0; i < fft_size; i++) {
        int bin_idx = static_cast<int>(power[i] / bin_width);
        if (bin_idx >= bins) bin_idx = bins - 1;
        bin_counts[bin_idx]++;
    }

    // Normalize bin counts to probabilities
    double total = fft_size;
    for (int i = 0; i < bins; i++) {
        bin_counts[i] /= total;
    }

    // Calculate Shannon entropy
    double entropy = 0.0;
    for (int i = 0; i < bins; i++) {
        if (bin_counts[i] > 0) {
            entropy -= bin_counts[i] * std::log(bin_counts[i]);
        }
    }

    return entropy;
}

// ============================================================================
// Supplementary tsfresh features (simple aggregations, duplicates, booleans)
// ============================================================================

// [[Rcpp::export]]
double cpp_abs_energy(NumericVector x) {
    /**
     * Absolute energy of time series
     *
     * Calculates: sum(x^2)
     *
     * @param x Numeric vector (time series)
     * @return Sum of squared values
     */
    int n = x.size();
    if (n == 0) return NA_REAL;

    double sum_sq = 0.0;
    int count = 0;

    for (int i = 0; i < n; i++) {
        if (!NumericVector::is_na(x[i])) {
            sum_sq += x[i] * x[i];
            count++;
        }
    }

    if (count == 0) return NA_REAL;
    return sum_sq;
}

// [[Rcpp::export]]
double cpp_sum_values(NumericVector x) {
    /**
     * Sum of all values in time series
     *
     * @param x Numeric vector (time series)
     * @return Sum of all values
     */
    int n = x.size();
    if (n == 0) return NA_REAL;

    double sum = 0.0;
    int count = 0;

    for (int i = 0; i < n; i++) {
        if (!NumericVector::is_na(x[i])) {
            sum += x[i];
            count++;
        }
    }

    if (count == 0) return NA_REAL;
    return sum;
}

// [[Rcpp::export]]
bool cpp_has_duplicate(NumericVector x) {
    /**
     * Check if time series has duplicate values
     *
     * @param x Numeric vector (time series)
     * @return True if any value appears more than once
     */
    int n = x.size();
    if (n <= 1) return false;

    std::unordered_set<double> seen;

    for (int i = 0; i < n; i++) {
        if (!NumericVector::is_na(x[i])) {
            if (seen.find(x[i]) != seen.end()) {
                return true;  // Found duplicate
            }
            seen.insert(x[i]);
        }
    }

    return false;
}

// [[Rcpp::export]]
bool cpp_has_duplicate_max(NumericVector x) {
    /**
     * Check if maximum value appears more than once
     *
     * @param x Numeric vector (time series)
     * @return True if max appears multiple times
     */
    int n = x.size();
    if (n <= 1) return false;

    // Find maximum
    double max_val = -std::numeric_limits<double>::infinity();
    bool found_max = false;

    for (int i = 0; i < n; i++) {
        if (!NumericVector::is_na(x[i])) {
            if (!found_max || x[i] > max_val) {
                max_val = x[i];
                found_max = true;
            }
        }
    }

    if (!found_max) return false;

    // Count occurrences of max
    int count = 0;
    for (int i = 0; i < n; i++) {
        if (!NumericVector::is_na(x[i]) && x[i] == max_val) {
            count++;
            if (count > 1) return true;
        }
    }

    return false;
}

// [[Rcpp::export]]
bool cpp_has_duplicate_min(NumericVector x) {
    /**
     * Check if minimum value appears more than once
     *
     * @param x Numeric vector (time series)
     * @return True if min appears multiple times
     */
    int n = x.size();
    if (n <= 1) return false;

    // Find minimum
    double min_val = std::numeric_limits<double>::infinity();
    bool found_min = false;

    for (int i = 0; i < n; i++) {
        if (!NumericVector::is_na(x[i])) {
            if (!found_min || x[i] < min_val) {
                min_val = x[i];
                found_min = true;
            }
        }
    }

    if (!found_min) return false;

    // Count occurrences of min
    int count = 0;
    for (int i = 0; i < n; i++) {
        if (!NumericVector::is_na(x[i]) && x[i] == min_val) {
            count++;
            if (count > 1) return true;
        }
    }

    return false;
}

// [[Rcpp::export]]
double cpp_percentage_reoccurring(NumericVector x) {
    /**
     * Percentage of non-unique values to all values
     *
     * Returns (n_total - n_unique) / n_total
     *
     * @param x Numeric vector (time series)
     * @return Percentage of reoccurring values (0 to 1)
     */
    int n = x.size();
    if (n == 0) return NA_REAL;

    std::unordered_set<double> unique_values;
    int valid_count = 0;

    for (int i = 0; i < n; i++) {
        if (!NumericVector::is_na(x[i])) {
            unique_values.insert(x[i]);
            valid_count++;
        }
    }

    if (valid_count == 0) return NA_REAL;

    int n_unique = unique_values.size();
    return static_cast<double>(valid_count - n_unique) / valid_count;
}

// [[Rcpp::export]]
double cpp_sum_reoccurring(NumericVector x) {
    /**
     * Sum of values that appear more than once
     *
     * @param x Numeric vector (time series)
     * @return Sum of reoccurring values
     */
    int n = x.size();
    if (n == 0) return NA_REAL;

    // Count occurrences of each value
    std::map<double, int> counts;

    for (int i = 0; i < n; i++) {
        if (!NumericVector::is_na(x[i])) {
            counts[x[i]]++;
        }
    }

    // Sum values that appear more than once
    double sum = 0.0;
    for (const auto& pair : counts) {
        if (pair.second > 1) {
            sum += pair.first * pair.second;
        }
    }

    return sum;
}

// [[Rcpp::export]]
double cpp_ratio_unique_values(NumericVector x) {
    /**
     * Ratio of unique values to time series length
     *
     * Returns n_unique / n_total
     *
     * @param x Numeric vector (time series)
     * @return Ratio (0 to 1)
     */
    int n = x.size();
    if (n == 0) return NA_REAL;

    std::unordered_set<double> unique_values;
    int valid_count = 0;

    for (int i = 0; i < n; i++) {
        if (!NumericVector::is_na(x[i])) {
            unique_values.insert(x[i]);
            valid_count++;
        }
    }

    if (valid_count == 0) return NA_REAL;

    return static_cast<double>(unique_values.size()) / valid_count;
}

// [[Rcpp::export]]
bool cpp_symmetry_looking(NumericVector x, double r = 0.05) {
    /**
     * Check if distribution looks symmetric
     *
     * Returns True if |mean - median| < r * std
     *
     * @param x Numeric vector (time series)
     * @param r Threshold ratio (default 0.05)
     * @return True if symmetric looking
     */
    int n = x.size();
    if (n < 3) return false;

    // Calculate mean
    double sum = 0.0;
    int count = 0;
    for (int i = 0; i < n; i++) {
        if (!NumericVector::is_na(x[i])) {
            sum += x[i];
            count++;
        }
    }

    if (count < 3) return false;
    double mean = sum / count;

    // Calculate median
    std::vector<double> valid_values;
    for (int i = 0; i < n; i++) {
        if (!NumericVector::is_na(x[i])) {
            valid_values.push_back(x[i]);
        }
    }
    std::sort(valid_values.begin(), valid_values.end());

    double median;
    int mid = valid_values.size() / 2;
    if (valid_values.size() % 2 == 0) {
        median = (valid_values[mid - 1] + valid_values[mid]) / 2.0;
    } else {
        median = valid_values[mid];
    }

    // Calculate standard deviation
    double sum_sq = 0.0;
    for (const auto& val : valid_values) {
        double diff = val - mean;
        sum_sq += diff * diff;
    }
    double std = std::sqrt(sum_sq / count);

    if (std < 1e-10) return true;  // Constant series is symmetric

    // Check symmetry condition
    return std::abs(mean - median) < r * std;
}

// [[Rcpp::export]]
bool cpp_large_standard_deviation(NumericVector x, double r = 0.25) {
    /**
     * Check if standard deviation is large relative to range
     *
     * Returns True if std > r * (max - min)
     *
     * @param x Numeric vector (time series)
     * @param r Threshold ratio (default 0.25)
     * @return True if std is large
     */
    int n = x.size();
    if (n < 2) return false;

    double min_val = std::numeric_limits<double>::infinity();
    double max_val = -std::numeric_limits<double>::infinity();
    double sum = 0.0;
    int count = 0;

    for (int i = 0; i < n; i++) {
        if (!NumericVector::is_na(x[i])) {
            if (x[i] < min_val) min_val = x[i];
            if (x[i] > max_val) max_val = x[i];
            sum += x[i];
            count++;
        }
    }

    if (count < 2) return false;

    double mean = sum / count;
    double range = max_val - min_val;

    if (range < 1e-10) return false;  // Constant series

    // Calculate standard deviation
    double sum_sq = 0.0;
    for (int i = 0; i < n; i++) {
        if (!NumericVector::is_na(x[i])) {
            double diff = x[i] - mean;
            sum_sq += diff * diff;
        }
    }
    double std = std::sqrt(sum_sq / count);

    return std > r * range;
}

// [[Rcpp::export]]
bool cpp_variance_larger_than_std(NumericVector x) {
    /**
     * Check if variance is larger than standard deviation
     *
     * Returns True if var > std (mathematically, this means var > 1)
     *
     * @param x Numeric vector (time series)
     * @return True if variance > standard deviation
     */
    int n = x.size();
    if (n < 2) return false;

    double sum = 0.0;
    int count = 0;

    for (int i = 0; i < n; i++) {
        if (!NumericVector::is_na(x[i])) {
            sum += x[i];
            count++;
        }
    }

    if (count < 2) return false;

    double mean = sum / count;

    // Calculate variance
    double sum_sq = 0.0;
    for (int i = 0; i < n; i++) {
        if (!NumericVector::is_na(x[i])) {
            double diff = x[i] - mean;
            sum_sq += diff * diff;
        }
    }
    double variance = sum_sq / count;

    // Variance > std means variance > sqrt(variance), i.e., variance > 1
    return variance > 1.0;
}
