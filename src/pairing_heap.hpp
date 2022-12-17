/*
 * An implementation of a minimum pairing heap for use with
 * 3D djikstra. Can be easily generalized.
 *
 * Michael L. Fredman, Robert Sedgewick, Daniel D. Sleator, 
 * and Robert E. Tarjan
 * "The pairing heap: A new form of self-adjusting heap."
 * Algorithmica. Nov. 1986, Vol. 1, Iss. 1-4, pp. 111-129
 * doi: 10.1007/BF01840439
 *
 *
 * Author: William Silversmith
 * Affiliation: Seung Lab, Princeton University
 * Date: August 2018
 */

#ifndef PAIRING_HEAP_HPP
#define PAIRING_HEAP_HPP

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <vector>

namespace crackle {
namespace pairing_heap {

template <typename KEY, typename VALUE>
class PHNode {
typedef PHNode<KEY, VALUE> PHNode_t;

public:
  PHNode_t* left;
  PHNode_t* right;
  KEY key; 
  VALUE value;

  // pp. 114: "In order to make "decrease key" and "delete"
  // more efficient, we must store with each node a third pointer, 
  // to its parent in the binary tree."

  PHNode_t* parent; 

  PHNode() {
    left = NULL;
    right = NULL;
    parent = NULL;
    key = 0;
    value = 0;
  }

  PHNode(KEY k, VALUE val) {
    left = NULL;
    right = NULL;
    parent = NULL;
    key = k;
    value = val;
  }

  PHNode(PHNode_t* lt, PHNode_t* rt, PHNode_t* p, KEY k, VALUE val) {
    left = lt;
    right = rt;
    parent = p;
    key = k;
    value = val;
  }

  PHNode (const PHNode_t &p) {
    left = p.left;
    right = p.right;
    parent = p.parent;
    key = p.key;
    value = p.value;
  }

  ~PHNode () {}
};

// O(1)
template <typename KEY, typename VALUE>
PHNode<KEY,VALUE>* meld(
  PHNode<KEY,VALUE>* h1, PHNode<KEY,VALUE>* h2
) {
  if (h1->key <= h2->key) {
    h2->right = h1->left;
    h1->left = h2;
    h2->parent = h1;
    return h1;
  }
  
  h1->right = h2->left;
  h2->left = h1;
  h1->parent = h2;
  
  return h2;
}

// O(log n) amortized?
template <typename KEY, typename VALUE>
PHNode<KEY,VALUE>* popmin (PHNode<KEY,VALUE>* root) {
  PHNode<KEY,VALUE> *subtree = root->left;
  
  if (!subtree) {
    return NULL;
  }

  std::vector<PHNode<KEY,VALUE>*> forest;
  forest.reserve(16);

  while (subtree) {
    forest.push_back(subtree);
    subtree = subtree->right;
  }

  const uint64_t forest_size = forest.size();

  for (uint64_t i = 0; i < forest_size; i++) {
    forest[i]->parent = NULL;
    forest[i]->right = NULL;
  }

  if (forest_size == 1) {
    return forest[0];
  }

  // need to deal with lone subtrees?

  // forward pass
  uint64_t last = forest_size & 0xfffffffe; // if odd, size - 1
  for (uint64_t i = 0; i < last; i += 2) {
    forest[i >> 1] = meld(forest[i], forest[i + 1]); 
  }
  last >>= 1;

  if (forest_size & 0x1) { // if odd
    forest[last] = forest[forest_size - 1];
  }
  else {
    last--;
  }

  // backward pass
  for (uint64_t i = last; i > 0; i--) {
    forest[i-1] = meld(forest[i], forest[i - 1]);
  }

  return forest[0];
}

template <typename KEY, typename VALUE>
PHNode<KEY,VALUE>* delmin (PHNode<KEY,VALUE>* root) {
  PHNode<KEY,VALUE>* tmp = popmin(root);
  delete root;
  return tmp;
}

template <typename KEY, typename VALUE>
void unlink_parent(PHNode<KEY,VALUE>* node) {
    if (node->parent == NULL) {
      return;
    }
    
    if (node->parent->left == node) {
      node->parent->left = node->right;
    }
    else {
      PHNode<KEY,VALUE>* sib = node->parent->left;
      while (sib->right != node) {
        sib = sib->right;
      }
      sib->right = node->right;
    }
    node->parent = NULL;
    node->right = NULL;
}

template <typename KEY, typename VALUE>
class MinHeap {
typedef PHNode<KEY, VALUE> PHNode_t;
private:
  uint64_t _size;
public:
  PHNode_t* root;

  MinHeap() {
    root = NULL;
    _size = 0;
  }

  MinHeap (KEY key, const VALUE val) {
    root = new PHNode_t(key, val);
    _size = 1;
  }

  // O(n)
  ~MinHeap() {
    recursive_delete(root);
    root = NULL;
  }

  uint64_t size() {
    return _size;
  }

  bool empty () {
    return root == NULL;
  }

  float min_key () {
    if (root) {
      return root->key;
    }

    throw std::runtime_error("No min key.");
  }

  uint32_t min_value () {
    if (root) {
      return root->value;
    }

    throw std::runtime_error("No min value.");
  }

  // O(1)
  PHNode_t* min () {
    return root;
  }

  // O(1)
  PHNode_t* insert(KEY key, const uint32_t val) {
    PHNode_t* I = new PHNode_t(key, val);
    return insert(I);
  }

  // O(1)
  PHNode_t* insert(PHNode_t* I) {
    if (I == NULL) {
      return I;
    }

    _size++;

    if (root == NULL) {
      root = I;
      return I;
    }

    if (root->key <= I->key) {
      I->right = root->left;
      root->left = I;
      I->parent = root;
    }
    else {
      root->right = I->left;
      I->left = root;
      root->parent = I;
      root = I;
    }

    return I;
  }

  // O(1)
  void update_key(PHNode_t* x, KEY key) {
    KEY oldkey = x->key;
    x->key = key;

    if (x == root || oldkey == key) {
      return;
    }
    else if (x->parent->key < key) {
      return;
    }
    
    unlink_parent<KEY,VALUE>(x);

    if (oldkey < key) {
      PHNode_t* subtree = popmin(x);
      x->left = NULL;
      x->right = NULL;

      if (subtree == NULL) {
        root = meld(x, root);
      }
      else {
        root = meld(x, meld(root, subtree));
      }
    }
    else {
      root = meld(x, root);
    }
  }

  // O(log n) amortized?
  void pop () {
    if (!root) {
      return;
    }
    _size--;
    
    root = delmin(root);
  }

  void erase (PHNode_t* x) {
    _size--;

    if (x == root) {
      root = delmin(root);
      return;
    } 

    unlink_parent<KEY,VALUE>(x);

    x = delmin(x);
    if (x != NULL) {
      root = meld(root, x);
    }
  }

private:
  void recursive_delete (PHNode_t* n) {
    if (n == NULL) {
      return;
    }

    if (n->left != NULL) {
      recursive_delete(n->left);
    }

    if (n->right != NULL) {
      recursive_delete(n->right);
    }

    if (n->parent) {
      if (n->parent->left == n) {
        n->parent->left = NULL;
      }
      else if (n->parent->right == n) {
        n->parent->right = NULL;
      }
    }

    delete n;
  }
};

};
};

#endif