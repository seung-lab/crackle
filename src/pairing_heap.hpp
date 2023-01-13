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
  PHNode_t* child;
  PHNode_t* next_sib;
  PHNode_t* prev_sib;
  KEY key; 
  VALUE value;

  // pp. 114: "In order to make "decrease key" and "delete"
  // more efficient, we must store with each node a third pointer, 
  // to its parent in the binary tree."

  PHNode_t* parent; 

  PHNode() {
    child = NULL;
    next_sib = NULL;
    prev_sib = NULL;
    parent = NULL;
    key = 0;
    value = 0;
  }

  PHNode(KEY k, VALUE val) {
    child = NULL;
    next_sib = NULL;
    prev_sib = NULL;
    parent = NULL;
    key = k;
    value = val;
  }

  PHNode (const PHNode_t &p) {
    child = p.child;
    next_sib = p.next_sib;
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
    h2->next_sib = h1->child;
    if (h2->next_sib) {
      h2->next_sib->prev_sib = h2;
    }
    h1->child = h2;
    h2->parent = h1;
    return h1;
  }
  
  h1->next_sib = h2->child;
  if (h1->next_sib) {
    h1->next_sib->prev_sib = h1;
  }
  h2->child = h1;
  h1->parent = h2;
  
  return h2;
}

// O(log n) amortized?
template <typename KEY, typename VALUE>
PHNode<KEY,VALUE>* popmin (PHNode<KEY,VALUE>* root) {
  PHNode<KEY,VALUE> *subtree = root->child;
  
  if (!subtree) {
    return NULL;
  }

  std::vector<PHNode<KEY,VALUE>*> forest;
  forest.reserve(16);

  while (subtree) {
    forest.push_back(subtree);
    subtree = subtree->next_sib;
  }

  const uint64_t forest_size = forest.size();

  for (uint64_t i = 0; i < forest_size; i++) {
    forest[i]->parent = NULL;
    forest[i]->next_sib = NULL;
    forest[i]->prev_sib = NULL;
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
    
    if (node->parent->child == node) {
      node->parent->child = node->next_sib;
      if (node->next_sib) {
        node->next_sib->prev_sib = NULL;
      }
    }
    else {
      PHNode<KEY,VALUE>* sib = node->prev_sib;
      node->prev_sib = node->next_sib;
      sib->next_sib = node->next_sib;
      if (node->next_sib) {
        node->next_sib->prev_sib = sib;
      }
    }
    node->parent = NULL;
    node->next_sib = NULL;
    node->prev_sib = NULL;
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

  KEY min_key () {
    if (root) {
      return root->key;
    }

    throw std::runtime_error("No min key.");
  }

  VALUE min_value () {
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
  PHNode_t* emplace(KEY key, const uint32_t val) {
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

    root = meld(root, I);
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
      x->child = NULL;
      x->next_sib = NULL;
      x->prev_sib = NULL;

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

    if (n->child != NULL) {
      recursive_delete(n->child);
    }

    if (n->next_sib != NULL) {
      recursive_delete(n->next_sib);
    }

    if (n->parent) {
      if (n->parent->child == n) {
        n->parent->child = NULL;
      }
      else if (n->parent->next_sib == n) {
        n->parent->next_sib = NULL;
      }
    }

    delete n;
  }
};

};
};

#endif