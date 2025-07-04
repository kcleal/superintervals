#!/bin/bash

# SuperIntervals R Package Build Script

set -e  # Exit on any error

echo "=== Cleaning old build files ==="
rm -f src/*.o src/*.so

echo "=== Generating Rcpp exports ==="
Rscript -e "library(Rcpp); compileAttributes('.')"

echo "=== Updating documentation ==="
Rscript -e "library(roxygen2); roxygenise()"

echo "=== Building package ==="
R CMD build .

echo "=== Running R CMD check ==="
TARBALL=$(ls -t superintervals_*.tar.gz | head -n1)
R CMD check "$TARBALL" --no-manual

# echo "=== Build complete ==="
# echo "Package tarball: $TARBALL"

echo "=== Installing ==="
R CMD INSTALL .