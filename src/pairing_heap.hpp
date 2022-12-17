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
#include <stdio.h>
#include <vector>
#include <unistd.h>

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

  void print () {
    printf("PHNode[%p](%.1f, %d, %p, %p, %p)\n", this, key, value, left, right, parent);
  }
};

template <typename KEY, typename VALUE>
void really_print_keys(PHNode<KEY,VALUE>* n, const int depth) {
  printf("(%d) %1.f \n", depth, n->key);

  if (depth > 20) {
    return;
  }

  if (n->left != NULL) {
    printf("L");
    really_print_keys(n->left, depth+1);
  }

  if (n->right != NULL) {
    printf("R");
    really_print_keys(n->right, depth+1);
  }    
}

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
PHNode<KEY,VALUE>* delmin (PHNode<KEY,VALUE>* root) {
  PHNode<KEY,VALUE> *subtree = root->left;
  
  if (!subtree) {
    delete root;
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
    delete root;
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

  delete root;
  return forest[0];
}


template <typename KEY, typename VALUE>
class MinHeap {
typedef PHNode<KEY, VALUE> PHNode_t;
public:
  PHNode_t* root;

  MinHeap() {
    root = NULL;
  }

  MinHeap (KEY key, const VALUE val) {
    root = new PHNode_t(key, val);
  }

  // // O(n)
  // ~MinHeap() {
  //   recursive_delete(root);
  //   root = NULL;
  // }

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
  PHNode_t* find_min () {
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

  void decrease_key (KEY delta) {
    if (root) {
      root->key -= delta;
    }
  }

  // O(1)
  void update_key (PHNode_t* x, KEY key) {
    x->key = key;

    if (x == root) {
      return;
    }
    
    // Assuming I do this right,
    // x not being root should be sufficent
    // to mean it has a parent.
    // if (x->parent) {

    if (x->parent->left == x) {
      x->parent->left = NULL;
    }
    else {
      x->parent->right = NULL;
    }

    insert(x);
  }

  // O(log n) amortized?
  void delete_min () {
    if (!root) {
      return;
    }

    root = delmin(root);
  }

  void delete_node (PHNode_t* x) {
    if (x == root) {
      root = delmin(root);
      return;
    } 

    if (x->parent->left == x) {
      x->parent->left = NULL;
    }
    else {
      x->parent->right = NULL;
    }

    // probably unnecessary line
    x->parent = NULL;

    insert(delmin(x));
  }

  void print_keys () {
    if (root) {
      really_print_keys(root, 0);
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