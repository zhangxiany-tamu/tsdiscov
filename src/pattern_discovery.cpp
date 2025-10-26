#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <limits>
using namespace Rcpp;

// Z-normalize a subsequence
std::vector<double> z_normalize(const std::vector<double>& x) {
    int n = x.size();
    if (n == 0) return x;

    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += x[i];
    mean /= n;

    double std = 0.0;
    for (int i = 0; i < n; i++) {
        double diff = x[i] - mean;
        std += diff * diff;
    }
    std = std::sqrt(std / n);

    if (std < 1e-10) {
        // Constant sequence - return zeros
        return std::vector<double>(n, 0.0);
    }

    std::vector<double> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = (x[i] - mean) / std;
    }

    return result;
}

// Compute Euclidean distance between two z-normalized subsequences
double euclidean_distance(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) return std::numeric_limits<double>::infinity();

    double sum = 0.0;
    for (size_t i = 0; i < a.size(); i++) {
        double diff = a[i] - b[i];
        sum += diff * diff;
    }

    return std::sqrt(sum);
}

// Simplified Matrix Profile computation
// Returns min, mean, max of matrix profile, plus motif and discord locations
// [[Rcpp::export]]
List cpp_matrix_profile(NumericVector x, int window_size = -1) {
    int n = x.size();

    // Default result
    List default_result = List::create(
        Named("matrix_profile_min") = NA_REAL,
        Named("matrix_profile_mean") = NA_REAL,
        Named("matrix_profile_max") = NA_REAL,
        Named("motif_score") = NA_REAL,
        Named("discord_score") = NA_REAL
    );

    // Need sufficient data
    if (n < 20) return default_result;

    // Default window size: 10% of series length, min 4, max 100
    if (window_size < 0) {
        window_size = std::max(4, std::min(100, n / 10));
    }

    // Ensure window size is valid
    if (window_size < 4 || window_size > n / 2) {
        return default_result;
    }

    // Convert to std::vector for easier manipulation
    std::vector<double> x_vec(n);
    for (int i = 0; i < n; i++) {
        if (NumericVector::is_na(x[i])) {
            return default_result;  // Can't handle NAs in matrix profile
        }
        x_vec[i] = x[i];
    }

    // Number of subsequences
    int n_subseq = n - window_size + 1;

    // Matrix profile: for each subsequence, find nearest neighbor distance
    std::vector<double> matrix_profile(n_subseq);

    // Compute matrix profile (simplified version - O(n^2))
    // For large series, this should use STOMP algorithm, but for feature extraction
    // we just need representative statistics

    // To speed up for long series, sample subsequences if needed
    int max_comparisons = 10000;  // Limit computations
    int sample_rate = 1;
    if (n_subseq * n_subseq > max_comparisons) {
        sample_rate = std::ceil(std::sqrt(static_cast<double>(n_subseq * n_subseq) / max_comparisons));
    }

    for (int i = 0; i < n_subseq; i += sample_rate) {
        // Extract and normalize subsequence i
        std::vector<double> subseq_i(window_size);
        for (int k = 0; k < window_size; k++) {
            subseq_i[k] = x_vec[i + k];
        }
        std::vector<double> norm_i = z_normalize(subseq_i);

        // Find nearest neighbor (excluding trivial matches)
        double min_dist = std::numeric_limits<double>::infinity();

        for (int j = 0; j < n_subseq; j += sample_rate) {
            // Skip self and trivial matches (overlapping subsequences)
            if (std::abs(i - j) < window_size / 4) continue;

            // Extract and normalize subsequence j
            std::vector<double> subseq_j(window_size);
            for (int k = 0; k < window_size; k++) {
                subseq_j[k] = x_vec[j + k];
            }
            std::vector<double> norm_j = z_normalize(subseq_j);

            // Compute distance
            double dist = euclidean_distance(norm_i, norm_j);
            if (dist < min_dist) {
                min_dist = dist;
            }
        }

        matrix_profile[i / sample_rate] = min_dist;
    }

    // Adjust size if we sampled
    matrix_profile.resize(std::min(static_cast<size_t>(n_subseq), matrix_profile.size()));

    // Compute statistics
    if (matrix_profile.empty()) return default_result;

    double mp_min = *std::min_element(matrix_profile.begin(), matrix_profile.end());
    double mp_max = *std::max_element(matrix_profile.begin(), matrix_profile.end());

    double mp_sum = 0.0;
    for (double val : matrix_profile) {
        mp_sum += val;
    }
    double mp_mean = mp_sum / matrix_profile.size();

    // Motif score: negative of minimum (lower MP = better motif)
    // Normalized by mean for interpretability
    double motif_score = mp_mean > 0 ? -mp_min / mp_mean : 0.0;

    // Discord score: maximum divided by mean (higher MP = anomaly)
    double discord_score = mp_mean > 0 ? mp_max / mp_mean : 0.0;

    return List::create(
        Named("matrix_profile_min") = mp_min,
        Named("matrix_profile_mean") = mp_mean,
        Named("matrix_profile_max") = mp_max,
        Named("motif_score") = motif_score,
        Named("discord_score") = discord_score
    );
}

// Friedrich Coefficients - Polynomial coefficients from Langevin model
// Models dx/dt = f(x) + noise, where f(x) is approximated as polynomial
// [[Rcpp::export]]
List cpp_friedrich_coefficients(NumericVector x, int max_order = 3) {
    int n = x.size();

    // Default result
    List default_result;
    for (int i = 0; i <= max_order; i++) {
        std::string name = "friedrich_coef_" + std::to_string(i);
        default_result[name] = NA_REAL;
    }

    // Need sufficient data
    if (n < 20) return default_result;

    // Remove NAs
    std::vector<double> x_clean;
    for (int i = 0; i < n; i++) {
        if (!NumericVector::is_na(x[i])) {
            x_clean.push_back(x[i]);
        }
    }

    if (x_clean.size() < 20) return default_result;
    n = x_clean.size();

    // Check for variation
    double mean_x = 0.0;
    for (int i = 0; i < n; i++) mean_x += x_clean[i];
    mean_x /= n;

    double var_x = 0.0;
    for (int i = 0; i < n; i++) {
        double diff = x_clean[i] - mean_x;
        var_x += diff * diff;
    }
    var_x /= (n - 1);

    if (var_x < 1e-10) return default_result;

    // Estimate drift: (x[t+1] - x[t]) / dt
    // Assume dt = 1 for simplicity
    std::vector<double> drift(n - 1);
    std::vector<double> x_mid(n - 1);

    for (int i = 0; i < n - 1; i++) {
        drift[i] = x_clean[i + 1] - x_clean[i];
        x_mid[i] = x_clean[i];  // Use x[t] as predictor
    }

    // Fit polynomial: drift ~ a₀ + a₁·x + a₂·x² + ... + aₚ·xᵖ
    // Using least squares regression

    int n_obs = n - 1;
    int n_coef = max_order + 1;  // Number of coefficients (including intercept)

    // Build design matrix X where X[i,j] = x_mid[i]^j
    std::vector<std::vector<double>> X(n_obs, std::vector<double>(n_coef));
    for (int i = 0; i < n_obs; i++) {
        double x_val = x_mid[i];
        double x_power = 1.0;
        for (int j = 0; j < n_coef; j++) {
            X[i][j] = x_power;
            x_power *= x_val;
        }
    }

    // Compute X'X
    std::vector<std::vector<double>> XtX(n_coef, std::vector<double>(n_coef, 0.0));
    for (int i = 0; i < n_coef; i++) {
        for (int j = 0; j < n_coef; j++) {
            double sum = 0.0;
            for (int k = 0; k < n_obs; k++) {
                sum += X[k][i] * X[k][j];
            }
            XtX[i][j] = sum;
        }
    }

    // Compute X'y
    std::vector<double> Xty(n_coef, 0.0);
    for (int i = 0; i < n_coef; i++) {
        double sum = 0.0;
        for (int k = 0; k < n_obs; k++) {
            sum += X[k][i] * drift[k];
        }
        Xty[i] = sum;
    }

    // Solve (X'X) β = X'y using Gaussian elimination
    // Simple implementation for small systems

    // Create augmented matrix [XtX | Xty]
    std::vector<std::vector<double>> aug(n_coef, std::vector<double>(n_coef + 1));
    for (int i = 0; i < n_coef; i++) {
        for (int j = 0; j < n_coef; j++) {
            aug[i][j] = XtX[i][j];
        }
        aug[i][n_coef] = Xty[i];
    }

    // Forward elimination
    for (int k = 0; k < n_coef; k++) {
        // Find pivot
        int max_row = k;
        double max_val = std::abs(aug[k][k]);
        for (int i = k + 1; i < n_coef; i++) {
            if (std::abs(aug[i][k]) > max_val) {
                max_val = std::abs(aug[i][k]);
                max_row = i;
            }
        }

        // Swap rows
        if (max_row != k) {
            std::swap(aug[k], aug[max_row]);
        }

        // Check for singularity
        if (std::abs(aug[k][k]) < 1e-10) {
            // Singular matrix - cannot solve
            return default_result;
        }

        // Eliminate below pivot
        for (int i = k + 1; i < n_coef; i++) {
            double factor = aug[i][k] / aug[k][k];
            for (int j = k; j <= n_coef; j++) {
                aug[i][j] -= factor * aug[k][j];
            }
        }
    }

    // Back substitution
    std::vector<double> beta(n_coef);
    for (int i = n_coef - 1; i >= 0; i--) {
        double sum = aug[i][n_coef];
        for (int j = i + 1; j < n_coef; j++) {
            sum -= aug[i][j] * beta[j];
        }
        beta[i] = sum / aug[i][i];
    }

    // Return coefficients
    List result;
    for (int i = 0; i <= max_order; i++) {
        std::string name = "friedrich_coef_" + std::to_string(i);
        result[name] = beta[i];
    }

    return result;
}

// Langevin Fixed Point - Find stable states where drift = 0
// Uses Newton-Raphson to find roots of polynomial f(x) = a0 + a1*x + a2*x^2 + a3*x^3
// [[Rcpp::export]]
List cpp_langevin_fixed_point(NumericVector x, int max_order = 3) {
    // First get Friedrich coefficients
    List friedrich = cpp_friedrich_coefficients(x, max_order);

    // Default result
    List default_result = List::create(
        Named("langevin_fixed_point") = NA_REAL,
        Named("langevin_max_fixed_point") = NA_REAL
    );

    // Check if Friedrich coefficients are valid
    bool valid = true;
    std::vector<double> coef(max_order + 1);
    for (int i = 0; i <= max_order; i++) {
        std::string name = "friedrich_coef_" + std::to_string(i);
        if (Rf_isNull(friedrich[name]) || NumericVector::is_na(as<double>(friedrich[name]))) {
            valid = false;
            break;
        }
        coef[i] = as<double>(friedrich[name]);
    }

    if (!valid) return default_result;

    // Evaluate polynomial and its derivative
    auto eval_poly = [&](double x_val) {
        double result = 0.0;
        double x_power = 1.0;
        for (int i = 0; i <= max_order; i++) {
            result += coef[i] * x_power;
            x_power *= x_val;
        }
        return result;
    };

    auto eval_poly_deriv = [&](double x_val) {
        double result = 0.0;
        double x_power = 1.0;
        for (int i = 1; i <= max_order; i++) {
            result += i * coef[i] * x_power;
            x_power *= x_val;
        }
        return result;
    };

    // Find fixed points using Newton-Raphson from multiple starting points
    // Try starting points based on data range
    double x_min = min(x);
    double x_max = max(x);
    double x_mean = mean(x);
    double x_range = x_max - x_min;

    std::vector<double> starting_points = {
        x_min, x_mean, x_max,
        x_min + 0.25 * x_range,
        x_min + 0.75 * x_range
    };

    std::vector<double> fixed_points;
    double tolerance = 1e-6;
    int max_iter = 50;

    for (double x0 : starting_points) {
        double x_curr = x0;
        bool converged = false;

        for (int iter = 0; iter < max_iter; iter++) {
            double f = eval_poly(x_curr);
            double df = eval_poly_deriv(x_curr);

            if (std::abs(df) < 1e-10) break;  // Derivative too small

            double x_next = x_curr - f / df;

            if (std::abs(x_next - x_curr) < tolerance) {
                converged = true;
                x_curr = x_next;
                break;
            }

            x_curr = x_next;

            // Keep within reasonable bounds
            if (x_curr < x_min - x_range || x_curr > x_max + x_range) {
                break;
            }
        }

        if (converged) {
            // Check if this is a new fixed point (not already found)
            bool is_new = true;
            for (double fp : fixed_points) {
                if (std::abs(x_curr - fp) < tolerance * 10) {
                    is_new = false;
                    break;
                }
            }

            if (is_new) {
                fixed_points.push_back(x_curr);
            }
        }
    }

    if (fixed_points.empty()) {
        return default_result;
    }

    // Return first fixed point and maximum fixed point
    double first_fp = fixed_points[0];
    double max_fp = *std::max_element(fixed_points.begin(), fixed_points.end());

    return List::create(
        Named("langevin_fixed_point") = first_fp,
        Named("langevin_max_fixed_point") = max_fp
    );
}

