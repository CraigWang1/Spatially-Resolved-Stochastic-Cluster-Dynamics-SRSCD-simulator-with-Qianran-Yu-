#ifndef UTIL_H
#define UTIL_H

#include "constants.h"
#include <vector>

using namespace std;

double volumeAtIndex(int);
double length(int);
double lengthf(int);
double lengthb(int);

template<typename T>
struct SegmentTree {
    int n_orig;
    int n; 
    std::vector<T> tree;
    long long op_count = 0;
    const long long RESYNC_THRESHOLD = 1000000;

    SegmentTree(int n_in) : n_orig(n_in) {
        n = 1;
        while (n < n_in) n <<= 1;
        tree.assign(2 * n, T(0));
    }

    // Wipes out rounding drift by rebuilding from leaves up
    void resync() {
        for (int i = n - 1; i > 0; --i) {
            tree[i] = tree[2 * i] + tree[2 * i + 1];
        }
        op_count = 0;
    }

    // Initial build logic
    void build_from_vector(const std::vector<T>& arr_in) {
        for (int i = 0; i < n_orig; ++i) tree[n + i] = arr_in[i];
        for (int i = n_orig; i < n; ++i) tree[n + i] = 0;
        resync();
    }

    // Optimized Iterative set_val with automatic resync
    void set_val(int idx, T val) {
        int tree_idx = idx + n;
        tree[tree_idx] = val;
        
        // Update the specific branch
        while (tree_idx > 1) {
            tree_idx >>= 1;
            tree[tree_idx] = tree[2 * tree_idx] + tree[2 * tree_idx + 1];
        }

        // Automatic resync to prevent numerical drift
        if (++op_count >= RESYNC_THRESHOLD) {
            resync();
        }
    }

    // Iterative Range Sum [ql..qr]
    T sum(int ql, int qr) {
        if (ql > qr) return T(0);
        T res = 0;
        int l = ql + n;
        int r = qr + n + 1; // Open interval [l, r)
        
        while (l < r) {
            if (l & 1) res += tree[l++];
            if (r & 1) res += tree[--r];
            l >>= 1;
            r >>= 1;
        }
        return res;
    }

    // Search for Gillespie: find smallest i such that sum(0..i) >= target
    int first_prefix_at_least(T target) {
        if (target <= 0) return 0;
        if (target >= tree[1]) return n_orig - 1; // In rare edge case this happens, return last index

        int node = 1;
        while (node < n) {
            node <<= 1;
            if (tree[node] < target) {
                target -= tree[node];
                node++;
            }
        }
        int res = node - n;
        return (res < n_orig) ? res : n_orig;
    }

    int first_prefix_at_least_from(int ql, T target) {
        T prefix_before = (ql > 0) ? sum(0, ql - 1) : T(0);
        return first_prefix_at_least(target + prefix_before);
    }

    // Helper to get raw values if your code expects an array/vector access
    T get_val(int idx) const {
        return tree[n + idx];
    }
};

#endif