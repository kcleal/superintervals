// File: bench_probe2d_new.cpp
// Compile: g++ -O3 -std=c++17 superintervals_nd.cpp -o bench2d
// Run: ./bench2d > probe2d_results.csv 2> probe2d_buildlog.txt
// Results summary
// -------------------------------
// Alg             Number of wins
// KDtree          44
// IntervalMap2D   21
// BVH             1
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <tuple>
#include <vector>

using namespace std;

static constexpr double WORLD_MIN = 0.0;
static constexpr double WORLD_MAX = 100.0;

template<typename F>
static inline double time_ms(F&& f){
    auto a = chrono::high_resolution_clock::now();
    f();
    auto b = chrono::high_resolution_clock::now();
    return chrono::duration<double, milli>(b - a).count();
}

static inline void make_square(double cx, double cy, double L,
                               double& x0, double& x1, double& y0, double& y1){
    double h = 0.5 * L;
    x0 = cx - h; x1 = cx + h;
    y0 = cy - h; y1 = cy + h;
    // clamp to [WORLD_MIN, WORLD_MAX]
    if (x0 < WORLD_MIN){ x1 -= (x0 - WORLD_MIN); x0 = WORLD_MIN; }
    if (x1 > WORLD_MAX){ x0 -= (x1 - WORLD_MAX); x1 = WORLD_MAX; }
    if (y0 < WORLD_MIN){ y1 -= (y0 - WORLD_MIN); y0 = WORLD_MIN; }
    if (y1 > WORLD_MAX){ y0 -= (y1 - WORLD_MAX); y1 = WORLD_MAX; }
    if (x0 < WORLD_MIN) x0 = WORLD_MIN;
    if (y0 < WORLD_MIN) y0 = WORLD_MIN;
}

// ============================================================================
// IntervalMap2D: superinterval traversal on X, linear filter on Y (SoA)
// Basically, this is a super-fast 1D interval index on X, and Y is just a cheap filter.
// This only works well if X alone already does most of the 'candidate culling' or if
// the query is big enough that candidate counts are high anyway.
// Use cases:
// - If selectivity is low (few hits per query), KD/BVH-like structures usually win.
// - If selectivity is moderate-to-high (many hits per query),
//   IntervalMapND becomes competitive and can win.
// - If fast build times are needed, IntervalMap2D is several times faster
// ============================================================================
class IntervalMap2D {
public:
    using S = double;
    vector<S> sx, ex, sy, ey;
    vector<size_t> branch;

    void reserve(size_t n){
        sx.reserve(n); ex.reserve(n); sy.reserve(n); ey.reserve(n);
        branch.reserve(n);
    }

    void add(S x0, S x1, S y0, S y1){
        sx.push_back(x0); ex.push_back(x1);
        sy.push_back(y0); ey.push_back(y1);
    }

    size_t size() const { return sx.size(); }

    void build(){
        size_t n = size();
        if(!n) return;

        vector<size_t> idx(n);
        iota(idx.begin(), idx.end(), 0);
        sort(idx.begin(), idx.end(), [&](size_t a, size_t b){
            if (sx[a] != sx[b]) return sx[a] < sx[b];
            return ex[a] > ex[b];
        });

        auto perm = [&](vector<S>& v){
            vector<S> out(n);
            for(size_t i=0;i<n;i++) out[i] = v[idx[i]];
            v.swap(out);
        };
        perm(sx); perm(ex); perm(sy); perm(ey);

        branch.assign(n, SIZE_MAX);
        vector<pair<S,size_t>> st;
        st.reserve(n);
        for(size_t i=0;i<n;i++){
            while(!st.empty() && st.back().first <= ex[i]) st.pop_back();
            if(!st.empty()) branch[i] = st.back().second;
            st.push_back({ex[i], i});
        }
    }

    inline size_t upper_bound_x(S value) const {
        size_t n = sx.size();
        if(!n) return SIZE_MAX;
        size_t len = n;
        size_t idx = 0;
        while(len > 1){
            size_t half = len / 2;
            idx += (sx[idx + half] <= value) * (len - half);
            len = half;
        }
        if (sx[idx] > value){
            if (idx == 0) return SIZE_MAX;
            --idx;
        }
        return idx;
    }

    size_t count(S qx0, S qx1, S qy0, S qy1) const {
        if (sx.empty()) return 0;
        size_t i = upper_bound_x(qx1);
        if (i == SIZE_MAX) return 0;
        size_t c = 0;
        while(i != SIZE_MAX){
            if (qx0 <= ex[i]){
                if (!(sy[i] > qy1 || qy0 > ey[i])) ++c;
                --i;
            } else {
                i = branch[i];
            }
        }
        return c;
    }
};

// ============================================================================
// KDTreeFast2D: flat node array, nth_element build, iterative query
// ============================================================================
class KDTreeFast2D {
public:
    using S = double;
    struct Node{
        int l=-1, r=-1;
        uint32_t bi=0;
        S mnx, mny, mxx, mxy;
        uint8_t split_dim=0;
    };

    vector<S> sx, ex, sy, ey;
    vector<uint32_t> idx;
    vector<Node> nodes;
    int root=-1;

    void reserve(size_t n){
        sx.reserve(n); ex.reserve(n); sy.reserve(n); ey.reserve(n);
        idx.reserve(n); nodes.reserve(n);
    }

    void add(S x0, S x1, S y0, S y1){
        sx.push_back(x0); ex.push_back(x1);
        sy.push_back(y0); ey.push_back(y1);
    }

    inline S center(uint32_t i, int d) const {
        return d==0 ? (sx[i]+ex[i])*0.5 : (sy[i]+ey[i])*0.5;
    }

    int build_rec(int l, int r, int depth){
        if (l>=r) return -1;
        int m = (l+r)/2;
        int d = depth & 1;

        nth_element(idx.begin()+l, idx.begin()+m, idx.begin()+r,
            [&](uint32_t a, uint32_t b){ return center(a,d) < center(b,d); });

        Node n;
        n.bi = idx[m];
        n.split_dim = (uint8_t)d;
        n.mnx = numeric_limits<S>::max();
        n.mny = numeric_limits<S>::max();
        n.mxx = numeric_limits<S>::lowest();
        n.mxy = numeric_limits<S>::lowest();

        for(int i=l;i<r;i++){
            uint32_t bi = idx[i];
            n.mnx = min(n.mnx, sx[bi]); n.mny = min(n.mny, sy[bi]);
            n.mxx = max(n.mxx, ex[bi]); n.mxy = max(n.mxy, ey[bi]);
        }

        int id = (int)nodes.size();
        nodes.push_back(n);
        nodes[id].l = build_rec(l, m, depth+1);
        nodes[id].r = build_rec(m+1, r, depth+1);
        return id;
    }

    void build(){
        size_t n = sx.size();
        idx.resize(n);
        iota(idx.begin(), idx.end(), 0);
        nodes.clear();
        nodes.reserve(n);
        root = build_rec(0, (int)n, 0);
    }

    size_t count(S qx0, S qx1, S qy0, S qy1) const {
        if (root < 0) return 0;
        size_t c=0;
        vector<int> st;
        st.reserve(64);
        st.push_back(root);
        while(!st.empty()){
            int ni = st.back(); st.pop_back();
            const Node& n = nodes[ni];
            if (n.mnx > qx1 || n.mxx < qx0 || n.mny > qy1 || n.mxy < qy0) continue;

            uint32_t bi = n.bi;
            if (!(sx[bi] > qx1 || qx0 > ex[bi] || sy[bi] > qy1 || qy0 > ey[bi])) ++c;
            if (n.l >= 0) st.push_back(n.l);
            if (n.r >= 0) st.push_back(n.r);
        }
        return c;
    }
};

// ============================================================================
// BVH2D: flat node array, median split along widest extent, iterative query
// ============================================================================
class BVH2D {
public:
    using S = double;
    struct Node{
        int l=-1, r=-1;
        uint32_t bi=0; // valid only for leaf
        S mnx, mny, mxx, mxy;
    };

    vector<S> sx, ex, sy, ey;
    vector<uint32_t> idx;
    vector<Node> nodes;
    int root=-1;

    void reserve(size_t n){
        sx.reserve(n); ex.reserve(n); sy.reserve(n); ey.reserve(n);
        idx.reserve(n); nodes.reserve(n*2);
    }

    void add(S x0, S x1, S y0, S y1){
        sx.push_back(x0); ex.push_back(x1);
        sy.push_back(y0); ey.push_back(y1);
    }

    int build_rec(int l, int r){
        if (l>=r) return -1;

        Node node;
        node.mnx = numeric_limits<S>::max();
        node.mny = numeric_limits<S>::max();
        node.mxx = numeric_limits<S>::lowest();
        node.mxy = numeric_limits<S>::lowest();

        for(int i=l;i<r;i++){
            uint32_t bi = idx[i];
            node.mnx = min(node.mnx, sx[bi]); node.mny = min(node.mny, sy[bi]);
            node.mxx = max(node.mxx, ex[bi]); node.mxy = max(node.mxy, ey[bi]);
        }

        int id = (int)nodes.size();
        nodes.push_back(node);

        if (r-l == 1){
            nodes[id].bi = idx[l];
            return id;
        }

        S exw = nodes[id].mxx - nodes[id].mnx;
        S eyw = nodes[id].mxy - nodes[id].mny;
        int dim = (eyw > exw) ? 1 : 0;
        int m = (l+r)/2;

        nth_element(idx.begin()+l, idx.begin()+m, idx.begin()+r, [&](uint32_t a, uint32_t b){
            double ca = dim==0 ? (sx[a]+ex[a]) : (sy[a]+ey[a]);
            double cb = dim==0 ? (sx[b]+ex[b]) : (sy[b]+ey[b]);
            return ca < cb;
        });

        int L = build_rec(l, m);
        int R = build_rec(m, r);
        nodes[id].l = L;
        nodes[id].r = R;
        return id;
    }

    void build(){
        size_t n = sx.size();
        idx.resize(n);
        iota(idx.begin(), idx.end(), 0);
        nodes.clear();
        nodes.reserve(n*2);
        root = (n ? build_rec(0, (int)n) : -1);
    }

    size_t count(S qx0, S qx1, S qy0, S qy1) const {
        if (root < 0) return 0;
        size_t c=0;
        vector<int> st;
        st.reserve(64);
        st.push_back(root);
        while(!st.empty()){
            int ni = st.back(); st.pop_back();
            const Node& n = nodes[ni];
            if (n.mnx > qx1 || n.mxx < qx0 || n.mny > qy1 || n.mxy < qy0) continue;

            if (n.l < 0 && n.r < 0){
                uint32_t bi = n.bi;
                if (!(sx[bi] > qx1 || qx0 > ex[bi] || sy[bi] > qy1 || qy0 > ey[bi])) ++c;
            } else {
                if (n.l >= 0) st.push_back(n.l);
                if (n.r >= 0) st.push_back(n.r);
            }
        }
        return c;
    }
};

int main(){
    const size_t n_boxes = 30000;
    const size_t n_queries = 500;
    const unsigned seed = 42;

    const vector<double> Lb_list = {2,5,10,20,30,40,60,80};
    const vector<double> Lq_list = {1,2,5,10,15,20,30,50};

    mt19937 rng(seed);
    uniform_real_distribution<double> dist(WORLD_MIN, WORLD_MAX);

    cout << "2D probe benchmark (squares)\n";
    cout << "n_boxes=" << n_boxes << ", n_queries=" << n_queries << ", world=[0,100]\n";
    cout << "Columns: Lb,Lq,selectivity_pct,IM_ms,KD_ms,BVH_ms,winner\\n\\n";
    cout << "Lb,Lq,selectivity_pct,IM_ms,KD_ms,BVH_ms,winner\\n";
    cout.setf(std::ios::fixed); cout<<setprecision(3);

    for (double Lb : Lb_list){
        IntervalMap2D im;
        KDTreeFast2D kd;
        BVH2D bvh;
        im.reserve(n_boxes);
        kd.reserve(n_boxes);
        bvh.reserve(n_boxes);

        for (size_t i=0;i<n_boxes;i++){
            double cx = dist(rng);
            double cy = dist(rng);
            double x0,x1,y0,y1;
            make_square(cx, cy, Lb, x0,x1,y0,y1);
            im.add(x0,x1,y0,y1);
            kd.add(x0,x1,y0,y1);
            bvh.add(x0,x1,y0,y1);
        }

        double im_build = time_ms([&]{ im.build(); });
        double kd_build = time_ms([&]{ kd.build(); });
        double bvh_build = time_ms([&]{ bvh.build(); });
        cerr << "Built Lb="<<Lb<<"  IM_build_ms="<<im_build
             << " KD_build_ms="<<kd_build<<" BVH_build_ms="<<bvh_build<<"\n";

        vector<pair<double,double>> qcenters;
        qcenters.reserve(n_queries);
        for (size_t i=0;i<n_queries;i++) qcenters.push_back({dist(rng), dist(rng)});

        for (double Lq : Lq_list){
            uint64_t hits_im=0, hits_kd=0, hits_bvh=0;

            vector<array<double,4>> queries;
            queries.reserve(n_queries);
            for (size_t i=0;i<n_queries;i++){
                double x0,x1,y0,y1;
                make_square(qcenters[i].first, qcenters[i].second, Lq, x0,x1,y0,y1);
                queries.push_back({x0,x1,y0,y1});
            }

            double t_im = time_ms([&]{
                for (size_t i=0;i<n_queries;i++){
                    const auto& q = queries[i];
                    hits_im += im.count(q[0],q[1],q[2],q[3]);
                }
            });

            double t_kd = time_ms([&]{
                for (size_t i=0;i<n_queries;i++){
                    const auto& q = queries[i];
                    hits_kd += kd.count(q[0],q[1],q[2],q[3]);
                }
            });

            double t_bvh = time_ms([&]{
                for (size_t i=0;i<n_queries;i++){
                    const auto& q = queries[i];
                    hits_bvh += bvh.count(q[0],q[1],q[2],q[3]);
                }
            });

            if (hits_im != hits_kd || hits_im != hits_bvh){
                cerr << "WARNING hit mismatch Lb="<<Lb<<" Lq="<<Lq
                     << " IM="<<hits_im<<" KD="<<hits_kd<<" BVH="<<hits_bvh<<"\n";
            }

            double select = 100.0 * (double)hits_im / ((double)n_boxes * (double)n_queries);

            string winner = "IM";
            double best = t_im;
            if (t_kd < best){ best = t_kd; winner = "KD"; }
            if (t_bvh < best){ best = t_bvh; winner = "BVH"; }

            cout << Lb << "," << Lq << "," << select << ","
                 << t_im << "," << t_kd << "," << t_bvh << ","
                 << winner << "\n";
        }
    }

    return 0;
}
