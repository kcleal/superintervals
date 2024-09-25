#include "superintervals.hpp"

template<typename S, typename T>
class SuperIntervalsE : public SuperIntervals<S, T> {
public:

    alignas(alignof(std::vector<S>)) std::vector<S> extent;

    void index() override {
        if (this->starts.size() == 0) {
            return;
        }
        this->starts.shrink_to_fit();
        this->ends.shrink_to_fit();
        this->data.shrink_to_fit();
        this->sortIntervals();

        eytz.resize(this->starts.size() + 1);
        eytz_index.resize(this->starts.size() + 1);
        eytzinger(&this->starts[0], this->starts.size());

        this->branch.resize(this->starts.size(), SIZE_MAX);
        std::vector<std::pair<S, size_t>> br;
        br.reserve(1000);
        br.emplace_back() = {this->ends[0], 0};
        for (size_t i=1; i < this->ends.size(); ++i) {
            while (!br.empty() && br.back().first < this->ends[i]) {
                br.pop_back();
            }
            if (!br.empty()) {
                this->branch[i] = br.back().second;
            }
            br.emplace_back() = {this->ends[i], i};
        }
        this->idx = 0;
    }

    inline void upperBound(const S x) noexcept override {
         size_t i = 0;
         const size_t n_intervals = this->starts.size();
         size_t best_idx = n_intervals;
         while (i < n_intervals) {
             if (eytz[i] > x) {
                 if (best_idx == n_intervals || eytz[i] <= eytz[best_idx]) {
                     best_idx = i;  // best candidate closer to x
                 }
                 i = 2 * i + 1;
             } else {
                 i = 2 * i + 2;
             }
         }
         // best_idx can be calculated afterwards like this, but turns out to be slower, not sure why:
         // int shift = __builtin_ffs(~(i + 1));
         // size_t best_idx = (i >> shift) - ((shift > 1) ? 1 : 0);
         this->idx = (best_idx < n_intervals) ? eytz_index[best_idx] : n_intervals - 1;
         if (this->idx > 0 && this->starts[this->idx] > x) {
             --this->idx;
         }
    }

private:
    std::vector<S> eytz;
    std::vector<size_t> eytz_index;

    size_t eytzinger_helper(S* arr, size_t n, size_t i, size_t k) {
        if (k < n) {
            i = eytzinger_helper(arr, n, i, 2*k+1);
            eytz[k] = this->starts[i];
            eytz_index[k] = i;
            ++i;
            i = eytzinger_helper(arr, n, i, 2*k + 2);
        }
        return i;
    }

    int eytzinger(S* arr, size_t n) {
        return eytzinger_helper(arr, n, 0, 0);
    }
};
