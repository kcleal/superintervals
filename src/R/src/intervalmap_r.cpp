#include <Rcpp.h>
#include "superintervals_r.h"

using namespace Rcpp;
using namespace si;

// Helper functions for SEXP reference management
inline void preserve_sexp(SEXP value) {
    if (value != R_NilValue) {
        R_PreserveObject(value);
    }
}

inline void release_sexp(SEXP value) {
    if (value != R_NilValue) {
        R_ReleaseObject(value);
    }
}

// Custom deleter for XPtr that handles SEXP cleanup
void cleanup_intervalmap(IntervalMap<int, SEXP>* ptr) {
    if (ptr) {
        // Release all SEXP references before destruction
        for (const SEXP& value : ptr->data) {
            release_sexp(value);
        }
        delete ptr;
    }
}

// [[Rcpp::export]]
SEXP create_intervalmap() {
    IntervalMap<int, SEXP>* ptr = new IntervalMap<int, SEXP>();
    XPtr<IntervalMap<int, SEXP>> main_ptr(ptr, &cleanup_intervalmap);

    // Create a cached search vectors and attach it to the XPtr as attributes
    std::vector<SEXP>* value_cache = new std::vector<SEXP>();
    XPtr<std::vector<SEXP>> value_cache_ptr(value_cache, true);
    main_ptr.attr("value_cache") = value_cache_ptr;

    std::vector<size_t>* idx_cache = new std::vector<size_t>();
    XPtr<std::vector<size_t>> idx_cache_ptr(idx_cache, true);
    main_ptr.attr("idx_cache") = idx_cache_ptr;

    return main_ptr;
}

// [[Rcpp::export]]
void add_interval(SEXP container, int start, int end, SEXP value = R_NilValue) {
    XPtr<IntervalMap<int, SEXP>> si(container);

    // Preserve the SEXP reference
    preserve_sexp(value);

    // Add to the interval map - this is now a simple operation
    si->add(start, end, value);
}

// [[Rcpp::export]]
void build_index(SEXP container) {
    XPtr<IntervalMap<int, SEXP>> si(container);
    si->build();
}

// [[Rcpp::export]]
List get_interval_at(SEXP container, int r_index) {
    XPtr<IntervalMap<int, SEXP>> si(container);

    // Convert R's 1-based to C++'s 0-based
    int cpp_index = r_index - 1;

    if (si->size() == 0 || cpp_index < 0 || cpp_index >= (int)si->size()) {
        stop("Index out of range");
    }

    int start = si->starts[cpp_index];
    int end = si->ends[cpp_index];
    SEXP value = si->data[cpp_index];

    return List::create(
        Named("start") = start,
        Named("end") = end,
        Named("value") = value
    );
}

// [[Rcpp::export]]
int get_start_at(SEXP container, int r_index) {
    XPtr<IntervalMap<int, SEXP>> si(container);

    // Convert R's 1-based to C++'s 0-based
    int cpp_index = r_index - 1;

    if (si->size() == 0 || cpp_index < 0 || cpp_index >= (int)si->size()) {
        stop("Index out of range");
    }

    return si->starts[cpp_index];
}

// [[Rcpp::export]]
int get_end_at(SEXP container, int r_index) {
    XPtr<IntervalMap<int, SEXP>> si(container);

    // Convert R's 1-based to C++'s 0-based
    int cpp_index = r_index - 1;

    if (si->size() == 0 || cpp_index < 0 || cpp_index >= (int)si->size()) {
        stop("Index out of range");
    }

    return si->ends[cpp_index];
}

// [[Rcpp::export]]
SEXP get_data_at(SEXP container, int r_index) {
    XPtr<IntervalMap<int, SEXP>> si(container);

    // Convert R's 1-based to C++'s 0-based
    int cpp_index = r_index - 1;

    if (si->size() == 0 || cpp_index < 0 || cpp_index >= (int)si->size()) {
        stop("Index out of range");
    }

    return si->data[cpp_index];
}

// [[Rcpp::export]]
void clear_intervals(SEXP container) {
    XPtr<IntervalMap<int, SEXP>> si(container);

    // Release all SEXP references before clearing
    for (const SEXP& value : si->data) {
        release_sexp(value);
    }

    si->clear();
}

// [[Rcpp::export]]
void reserve_intervals(SEXP container, int n) {
    XPtr<IntervalMap<int, SEXP>> si(container);
    si->reserve(n);
}

// [[Rcpp::export]]
int get_size(SEXP container) {
    XPtr<IntervalMap<int, SEXP>> si(container);
    return si->size();
}

// [[Rcpp::export]]
bool cpp_has_overlaps(SEXP container, int start, int end) {
    XPtr<IntervalMap<int, SEXP>> si(container);
    return si->has_overlaps(start, end);
}

// [[Rcpp::export]]
int count_overlaps(SEXP container, int start, int end) {
    XPtr<IntervalMap<int, SEXP>> si(container);
    return si->count(start, end);
}

// [[Rcpp::export]]
List cpp_search_values(SEXP container, int start, int end) {
    XPtr<IntervalMap<int, SEXP>> si(container);
    XPtr<std::vector<SEXP>> cache_ptr(si.attr("value_cache"));
    std::vector<SEXP>& found_values = *cache_ptr;
    found_values.clear();

    si->search_values(start, end, found_values);

    List result(found_values.size());
    for (size_t i = 0; i < found_values.size(); ++i) {
        result[i] = found_values[i];
    }

    return result;
}

// [[Rcpp::export]]
IntegerVector search_indexes(SEXP container, int start, int end) {
    XPtr<IntervalMap<int, SEXP>> si(container);
    XPtr<std::vector<size_t>> cache_ptr(si.attr("idx_cache"));
    std::vector<size_t>& found_indices = *cache_ptr;
    found_indices.clear();

    si->search_idxs(start, end, found_indices);

    IntegerVector result(found_indices.size());
    for (size_t i = 0; i < found_indices.size(); ++i) {
        result[i] = found_indices[i] + 1;  // Convert to 1-based indexing for R
    }

    return result;
}

// [[Rcpp::export]]
List cpp_search_keys(SEXP container, int start, int end) {
    XPtr<IntervalMap<int, SEXP>> si(container);
    XPtr<std::vector<size_t>> cache_ptr(si.attr("idx_cache"));
    std::vector<size_t>& found_indices = *cache_ptr;
    found_indices.clear();

    si->search_idxs(start, end, found_indices);

    if (found_indices.empty()) {
        return List::create(
            Named("start") = IntegerVector(),
            Named("end") = IntegerVector()
        );
    }

    IntegerVector starts_vec(found_indices.size());
    IntegerVector ends_vec(found_indices.size());

    for (size_t i = 0; i < found_indices.size(); ++i) {
        size_t idx = found_indices[i];
        starts_vec[i] = si->starts[idx];
        ends_vec[i] = si->ends[idx];
    }

    return List::create(
        Named("start") = starts_vec,
        Named("end") = ends_vec
    );
}

// [[Rcpp::export]]
List cpp_search_items(SEXP container, int start, int end) {
    XPtr<IntervalMap<int, SEXP>> si(container);
    XPtr<std::vector<size_t>> cache_ptr(si.attr("idx_cache"));
    std::vector<size_t>& found_indices = *cache_ptr;
    found_indices.clear();

    si->search_idxs(start, end, found_indices);

    if (found_indices.empty()) {
        return List::create(
            Named("start") = IntegerVector(),
            Named("end") = IntegerVector(),
            Named("value") = List()
        );
    }

    IntegerVector starts_vec(found_indices.size());
    IntegerVector ends_vec(found_indices.size());
    List values_list(found_indices.size());

    for (size_t i = 0; i < found_indices.size(); ++i) {
        size_t idx = found_indices[i];
        starts_vec[i] = si->starts[idx];
        ends_vec[i] = si->ends[idx];
        values_list[i] = si->data[idx];
    }

    return List::create(
        Named("start") = starts_vec,
        Named("end") = ends_vec,
        Named("value") = values_list
    );
}

// [[Rcpp::export]]
List get_coverage(SEXP container, int start, int end) {
    XPtr<IntervalMap<int, SEXP>> si(container);

    std::pair<size_t, int> cov_result;
    si->coverage(start, end, cov_result);

    return List::create(
        Named("count") = (int)cov_result.first,
        Named("total_coverage") = cov_result.second
    );
}