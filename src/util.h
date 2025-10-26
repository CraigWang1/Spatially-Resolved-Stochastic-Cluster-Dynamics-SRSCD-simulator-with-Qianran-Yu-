#ifndef UTIL_H
#define UTIL_H

#include "constants.h"
#include <unordered_map>
#include <vector>
#include <utility>

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

template <typename KeyT, typename ValueT>
class BinarySearchTree {
private:
    struct Node {
        KeyT id;
        ValueT rate;
        ValueT subtree_sum;
        Node* left;
        Node* right;

        Node(const KeyT& id_, const ValueT& rate_)
            : id(id_), rate(rate_), subtree_sum(rate_), left(nullptr), right(nullptr) {}
    };

    Node* root;
    bool dirty;
    std::vector<std::pair<KeyT, ValueT> > pending_inserts;
    std::vector<KeyT> pending_erases;
    std::unordered_map<KeyT, Node*> node_map; // O(1) lookup

    // --- Utility functions ---
    static ValueT get_sum(Node* node) {
        return node ? node->subtree_sum : ValueT();
    }

    static void clear(Node* node) {
        if (!node) return;
        clear(node->left);
        clear(node->right);
        delete node;
    }

    void recompute_all(Node* node) {
        if (!node) return;
        recompute_all(node->left);
        recompute_all(node->right);
        node->subtree_sum = get_sum(node->left) + node->rate + get_sum(node->right);
    }

    Node* insert(Node* node, const KeyT& id, const ValueT& rate) {
        if (!node) {
            Node* n = new Node(id, rate);
            node_map[id] = n;
            return n;
        }
        if (id < node->id) {
            node->left = insert(node->left, id, rate);
        } else if (id > node->id) {
            node->right = insert(node->right, id, rate);
        } else {
            node->rate = rate;
        }
        node->subtree_sum = get_sum(node->left) + node->rate + get_sum(node->right);
        return node;
    }

    Node* erase(Node* node, const KeyT& id, bool& erased) {
        if (!node) return nullptr;

        if (id < node->id) {
            node->left = erase(node->left, id, erased);
        } else if (id > node->id) {
            node->right = erase(node->right, id, erased);
        } else {
            erased = true;
            node_map.erase(id);
            if (!node->left) {
                Node* r = node->right;
                delete node;
                return r;
            } else if (!node->right) {
                Node* l = node->left;
                delete node;
                return l;
            } else {
                Node* succ = node->right;
                while (succ->left) succ = succ->left;
                node_map.erase(succ->id);
                node->id = succ->id;
                node->rate = succ->rate;
                node->right = erase(node->right, succ->id, erased = false);
                node_map[node->id] = node;
            }
        }
        node->subtree_sum = get_sum(node->left) + node->rate + get_sum(node->right);
        return node;
    }

    Node* find_node(const KeyT& id) const {
        typename std::unordered_map<KeyT, Node*>::const_iterator it = node_map.find(id);
        return it != node_map.end() ? it->second : nullptr;
    }

public:
    BinarySearchTree() : root(nullptr), dirty(false) {}

    ~BinarySearchTree() {
        clear(root);
    }

    // --- Batch update interface ---
    void batch_update_rate(const KeyT& id, const ValueT& new_rate) {
        Node* n = find_node(id);
        if (n) {
            n->rate = new_rate;
            dirty = true;
        } else {
            pending_inserts.push_back(std::make_pair(id, new_rate));
        }
    }

    void batch_erase(const KeyT& id) {
        if (dirty) {
            pending_erases.push_back(id);
        } else {
            bool erased = false;
            root = erase(root, id, erased);
        }
    }

    void insert_safe(const KeyT& id, const ValueT& rate) {
        if (dirty) {
            recompute_all(root);
            dirty = false;
        }

        // Apply any pending inserts or erases before inserting
        for (auto& e : pending_erases) {
            bool erased = false;
            root = erase(root, e, erased);
        }
        pending_erases.clear();

        for (auto& p : pending_inserts) {
            root = insert(root, p.first, p.second);
        }
        pending_inserts.clear();

        root = insert(root, id, rate);
    }

    void finalize_batch() {
        if (dirty) {
            recompute_all(root);
            dirty = false;
        }

        // apply pending erases
        for (typename std::vector<KeyT>::iterator it = pending_erases.begin(); it != pending_erases.end(); ++it) {
            bool erased = false;
            root = erase(root, *it, erased);
        }
        pending_erases.clear();

        // apply pending inserts
        for (typename std::vector<std::pair<KeyT, ValueT> >::iterator it = pending_inserts.begin(); it != pending_inserts.end(); ++it) {
            root = insert(root, it->first, it->second);
        }
        pending_inserts.clear();
    }

    KeyT find_prefix_ge(ValueT &sample) const {
        if (sample > get_total_rate()) {
            sample -= get_total_rate();
            return KeyT(0);
        }
        Node* node = root;
        while (node) {
            ValueT left_sum = get_sum(node->left);
            if (sample <= left_sum) {
                node = node->left;
            } else if (sample <= left_sum + node->rate) {
                sample -= left_sum;
                return node->id;
            } else {
                sample -= (left_sum + node->rate);
                node = node->right;
            }
        }

        // Shouldn't get here
        return KeyT(0);
    }

    ValueT get_total_rate() const {
        return get_sum(root);
    }
};


#endif