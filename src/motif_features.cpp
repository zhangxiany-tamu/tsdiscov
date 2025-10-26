#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
double cpp_motif_three_quantile(NumericVector x) {
    // SB_MotifThree_quantile_hh: symbolic analysis of 3-element patterns
    int n = x.size();
    if (n < 3) return NA_REAL;

    // Calculate quantiles for discretization
    NumericVector sorted = clone(x);
    std::sort(sorted.begin(), sorted.end());

    double q33 = sorted[n / 3];
    double q66 = sorted[2 * n / 3];

    // Convert to symbolic sequence (L=low, M=medium, H=high)
    std::vector<char> symbols(n);
    for (int i = 0; i < n; i++) {
        if (x[i] < q33) {
            symbols[i] = 'L';
        } else if (x[i] < q66) {
            symbols[i] = 'M';
        } else {
            symbols[i] = 'H';
        }
    }

    // Count 3-element patterns
    std::map<std::string, int> pattern_counts;
    int total_patterns = 0;

    for (int i = 0; i < n - 2; i++) {
        std::string pattern;
        pattern += symbols[i];
        pattern += symbols[i + 1];
        pattern += symbols[i + 2];
        pattern_counts[pattern]++;
        total_patterns++;
    }

    if (total_patterns == 0) return NA_REAL;

    // Calculate entropy of pattern distribution
    double entropy = 0.0;
    for (const auto& pair : pattern_counts) {
        double p = static_cast<double>(pair.second) / total_patterns;
        if (p > 0) {
            entropy -= p * std::log(p);
        }
    }

    // Focus on HH patterns (high-high transitions)
    std::map<std::string, int> hh_patterns;
    int hh_count = 0;

    for (int i = 0; i < n - 2; i++) {
        if (symbols[i] == 'H' && symbols[i + 1] == 'H') {
            std::string pattern;
            pattern += symbols[i];
            pattern += symbols[i + 1];
            pattern += symbols[i + 2];
            hh_patterns[pattern]++;
            hh_count++;
        }
    }

    if (hh_count == 0) return 0.0;

    // Calculate entropy of HH patterns
    double hh_entropy = 0.0;
    for (const auto& pair : hh_patterns) {
        double p = static_cast<double>(pair.second) / hh_count;
        if (p > 0) {
            hh_entropy -= p * std::log(p);
        }
    }

    return hh_entropy / std::log(2.0); // Normalize to bits
}
