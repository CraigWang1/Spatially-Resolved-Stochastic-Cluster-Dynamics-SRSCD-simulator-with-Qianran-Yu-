#ifndef UTIL_H
#define UTIL_H

#include "constants.h"
#include <vector>
#include <set>
#include <unordered_map>
#include <iostream>

using namespace std;

double volumeAtIndex(int);

template<typename T>
struct SegmentTree {
    int n;
    vector<T> tree;
    vector<T> arr;

    SegmentTree(int n) : n(n), tree(4*n, T()), arr(n, T()) {}

    SegmentTree() : SegmentTree(0) {}

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

    int size() {
    	return n;
    }

	void double_size() {
		if (n == 0)
			n = 1;
		else
			n *= 2;                       // double logical size
	    arr.resize(n, T());               // extend arr with default values
	    tree.assign(4*n, T());            // allocate new tree array
	    build(1, 0, n-1);                 // rebuild with new size
	}
};

/* Hash function for multisets */
struct MultisetHash {
    template <typename T>
    size_t operator()(const std::multiset<T>& ms) const {
        size_t h = 0;
        for (const auto& x : ms) {
            // boost-like hash combine
            h ^= std::hash<T>{}(x) + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        return h;
    }
};


template <typename K, typename V>
class BiMap {
    std::unordered_map<K, V> forward;
    std::unordered_map<V, K> backward;

public:
    bool insert(const K& k, const V& v) {
        if (forward.count(k) || backward.count(v)) return false; // enforce one-to-one
        forward[k] = v;
        backward[v] = k;
        return true;
    }

    const V& at_forward(const K& k) const { return forward.at(k); }
    const K& at_backward(const V& v) const { return backward.at(v); }

    void erase_by_key(const K& k) {
        auto it = forward.find(k);
        if (it != forward.end()) {
            backward.erase(it->second);
            forward.erase(it);
        }
    }

    void erase_by_value(const V& v) {
        auto it = backward.find(v);
        if (it != backward.end()) {
            forward.erase(it->second);
            backward.erase(it);
        }
    }

    bool contains_key(const K& k) const { return forward.count(k) > 0; }
    bool contains_value(const V& v) const { return backward.count(v) > 0; }
};


#endif