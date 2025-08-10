#ifndef UTIL_H
#define UTIL_H

#include "constants.h"
#include <vector>

using namespace std;

double volumeAtIndex(int);

template<typename T>
struct SegmentTree {
    int n;
    vector<T> tree;
    vector<T> arr;

    SegmentTree(int n) : n(n), tree(4*n, T()), arr(n, T()) {}

    void build(int node, int l, int r) {
        if (l == r) tree[node] = arr[l];
        else {
            int mid = (l + r) / 2;
            build(node*2, l, mid);
            build(node*2+1, mid+1, r);
            tree[node] = tree[node*2] + tree[node*2+1];
        }
    }

    // internal recursive set
    void set_val(int idx, T val, int node, int l, int r) {
        if (l == r) {
            arr[idx] = val;
            tree[node] = val;
            return;
        }
        int mid = (l + r) / 2;
        if (idx <= mid) set_val(idx, val, node*2, l, mid);
        else set_val(idx, val, node*2+1, mid+1, r);
        tree[node] = tree[node*2] + tree[node*2+1];
    }
    void set_val(int idx, T val) { set_val(idx, val, 1, 0, n-1); }

    // range sum [ql..qr]
    T sum(int ql, int qr, int node, int l, int r) {
        if (qr < l || ql > r) return T();
        if (ql <= l && r <= qr) return tree[node];
        int mid = (l + r) / 2;
        return sum(ql, qr, node*2, l, mid) + sum(ql, qr, node*2+1, mid+1, r);
    }
    T sum(int l, int r) {
        if (l > r) return T();
        return sum(l, r, 1, 0, n-1);
    }

    // --- search that assumes we're searching from index 0 ---
    // find smallest i such that sum(0..i) >= target
    int first_prefix_at_least_from_zero_rec(T target, int node, int l, int r) {
        if (l == r) return l;
        int mid = (l + r) / 2;
        if (tree[node*2] >= target) return first_prefix_at_least_from_zero_rec(target, node*2, l, mid);
        else return first_prefix_at_least_from_zero_rec(target - tree[node*2], node*2+1, mid+1, r);
    }
    // public wrapper (returns n if not found)
    int first_prefix_at_least(T target) {
        if (tree[1] < target) return n;
        return first_prefix_at_least_from_zero_rec(target, 1, 0, n-1);
    }

    // --- corrected: search starting from arbitrary ql ---
    // returns n when not found (you can change to -1 if you prefer)
    int first_prefix_at_least_from(int ql, T target) {
        if (ql < 0 || ql >= n) return n;
        if (target <= T()) { // if target <= 0, the first index is ql
            return ql;
        }
        T base = (ql == 0 ? T() : sum(0, ql-1));    // prefix sum up to ql-1
        T needed = base + target;                   // we need prefix >= needed
        if (tree[1] < needed) return n;             // overall total too small
        return first_prefix_at_least(needed);       // reuse the from-zero search
    }
};

#endif