# CRAN Submission Comments

## Test environments

* macOS Sequoia 15.6.1, R 4.5.0 (local)
* GitHub Actions (ubuntu-latest, windows-latest, macOS-latest), R-release

## R CMD check results

There were no ERRORs or WARNINGs.

There were 3 NOTEs:

1. **New submission**

   This is the first submission of tsdiscov to CRAN.

2. **pandoc not installed** (local check only)

   Local environment note: "Files 'README.md' or 'NEWS.md' cannot be checked
   without 'pandoc' being installed." CRAN's check servers have pandoc installed.

3. **HTML Tidy version** (local check only)

   Local environment note about HTML Tidy version. CRAN's check servers have
   the appropriate tools.

## Downstream dependencies

There are currently no downstream dependencies for this package.

## Additional information

* All 6003 tests pass without warnings or errors
* The package includes extensive C++ implementations using Rcpp and RcppArmadillo
* Comprehensive documentation is provided for all exported functions
* Examples are provided for all user-facing functions
