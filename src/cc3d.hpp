/*
 * Connected Components for 2D images. 

 * Author: William Silversmith
 * Affiliation: Seung Lab, Princeton University
 * Date: August 2018 - June 2019, 2021, Dec. 2022
 *
 * ----
 * LICENSE
 * 
 * This is a special reduced feature version of cc3d 
 * that includes only the logic needed for CCL 4-connected.
 * cc3d is ordinarily licensed as GPL v3. Get the full
 * version of cc3d here: 
 * 
 * https://github.com/seung-lab/connected-components-3d
 *
 *  License: GNU Lesser General Public License v3
 *  Copyright 2021 William Silversmith
 * 
 *  See COPYING.LESSER for that license. If the file is missing,
 *  you can find the text here: 
 *  https://www.gnu.org/licenses/lgpl-3.0.en.html
 */

#ifndef CC3D_SPECIAL_2_4_HPP
#define CC3D_SPECIAL_2_4_HPP 

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <stdexcept>

namespace cc3d {

static size_t _dummy_N;

template <typename T>
class DisjointSet {
public:
  T *ids;
  size_t length;

  DisjointSet () {
    length = 65536; // 2^16, some "reasonable" starting size
    ids = new T[length]();
  }

  DisjointSet (size_t len) {
    length = len;
    ids = new T[length]();
  }

  DisjointSet (const DisjointSet &cpy) {
    length = cpy.length;
    ids = new T[length]();

    for (int i = 0; i < length; i++) {
      ids[i] = cpy.ids[i];
    }
  }

  ~DisjointSet () {
    delete []ids;
  }

  T root (T n) {
    T i = ids[n];
    while (i != ids[i]) {
      ids[i] = ids[ids[i]]; // path compression
      i = ids[i];
    }

    return i;
  }

  bool find (T p, T q) {
    return root(p) == root(q);
  }

  void add(T p) {
    if (p >= length) {
      printf("Connected Components Error: Label %d cannot be mapped to union-find array of length %lu.\n", p, length);
      throw "maximum length exception";
    }

    if (ids[p] == 0) {
      ids[p] = p;
    }
  }

  void unify (T p, T q) {
    if (p == q) {
      return;
    }

    T i = root(p);
    T j = root(q);

    if (i == 0) {
      add(p);
      i = p;
    }

    if (j == 0) {
      add(q);
      j = q;
    }

    ids[i] = j;
  }
};

// This is the second raster pass of the two pass algorithm family.
// The input array (output_labels) has been assigned provisional 
// labels and this resolves them into their final labels. We
// modify this pass to also ensure that the output labels are
// numbered from 1 sequentially.
template <typename OUT = uint32_t>
OUT* relabel(
    OUT* out_labels, const int64_t voxels,
    const int64_t num_labels, DisjointSet<uint32_t> &equivalences,
    size_t &N = _dummy_N, OUT start_label = 1
  ) {

  OUT label;
  std::unique_ptr<OUT[]> renumber(new OUT[num_labels + 1]());
  OUT next_label = start_label;

  for (int64_t i = 1; i <= num_labels; i++) {
    label = equivalences.root(i);
    if (renumber[label] == 0) {
      renumber[label] = next_label;
      renumber[i] = next_label;
      next_label++;
    }
    else {
      renumber[i] = renumber[label];
    }
  }

  // Raster Scan 2: Write final labels based on equivalences
  N = next_label - start_label;
  if (N < static_cast<size_t>(num_labels) || start_label != 1) {
    for (int64_t loc = 0; loc < voxels; loc++) {
      out_labels[loc] = renumber[out_labels[loc]];
    }
  }

  return out_labels;
}

template <typename OUT>
OUT* color_connectivity_graph(
  uint8_t* vcg, // voxel connectivity graph
  const int64_t sx, const int64_t sy, const int64_t sz,
  size_t max_labels, OUT *out_labels = NULL, 
  size_t &N = _dummy_N, OUT start_label = 0
) {

  const int64_t sxy = sx * sy;
  const int64_t voxels = sx * sy * sz;

  max_labels++;
  max_labels = std::min(max_labels, static_cast<size_t>(voxels) + 1); // + 1L for an array with no zeros
  max_labels = std::min(max_labels, static_cast<size_t>(std::numeric_limits<OUT>::max()));

  DisjointSet<OUT> equivalences(max_labels);

  LABEL new_label = 0;
  equivalences.add(new_label);
  for (int64_t x = 0; x < sx; x++) {
    if (x > 0 && (vcg[x] & 0b0010) == 0) {
      new_label++;
      equivalences.add(new_label);
    }
    out[x] = new_label;
  }

  const int64_t B = -1;
  const int64_t C = -sx;
  const int64_t D = -1-sx;

  for (int64_t y = 1; y < sy; y++) {
    for (int64_t x = 0; x < sx; x++) {
      int64_t loc = x + sx * y;

      if (x > 0 && (vcg[loc] & 0b0010)) {
        out[loc] = out[loc+B];
        if (y > 0 && (vcg[loc + C] & 0b0010) == 0 && (vcg[loc] & 0b1000)) {
          equivalences.unify(out[loc], out[loc+C]);
        }
      }
      else if (y > 0 && vcg[loc] & 0b1000) {
        equivalences.unify(out[loc], out[loc+C]);
      }
      else {
        new_label++;
        out[loc] = new_label;
        equivalences.add(new_label);
      }
    }
  }

  return relabel<OUT>(out_labels, voxels, next_label, equivalences, N, start_label);
}

template <typename LABEL, typename OUT>
OUT* connected_components2d_4(
    LABEL* in_labels, 
    const int64_t sx, const int64_t sy, const int64_t sz,
    size_t max_labels, OUT *out_labels = NULL, 
    size_t &N = _dummy_N, OUT start_label = 0
  ) {

  const int64_t sxy = sx * sy;
  const int64_t voxels = sx * sy * sz;

  max_labels++;
  max_labels = std::min(max_labels, static_cast<size_t>(voxels) + 1); // + 1L for an array with no zeros
  max_labels = std::min(max_labels, static_cast<size_t>(std::numeric_limits<OUT>::max()));


  DisjointSet<OUT> equivalences(max_labels);

  if (out_labels == NULL) {
    out_labels = new OUT[voxels]();
  }
    
  /*
    Layout of forward pass mask. 
    A is the current location.
    D C 
    B A 
  */

  // const int64_t A = 0;
  const int64_t B = -1;
  const int64_t C = -sx;
  const int64_t D = -1-sx;

  int64_t loc = 0;
  OUT next_label = 0;

  LABEL cur = 0;
  for (int64_t z = 0; z < sz; z++) {
    for (int64_t y = 0; y < sy; y++) {
      for (int64_t x = 0; x < sx; x++) {
        loc = x + sx * y + sxy * z;
        cur = in_labels[loc];

        if (x > 0 && cur == in_labels[loc + B]) {
          out_labels[loc] = out_labels[loc + B];
          if (y > 0 && cur == in_labels[loc + C] && cur != in_labels[loc + D]) {
            equivalences.unify(out_labels[loc], out_labels[loc + C]);
          }
        }
        else if (y > 0 && cur == in_labels[loc + C]) {
          out_labels[loc] = out_labels[loc + C];
        }
        else {
          next_label++;
          out_labels[loc] = next_label;
          equivalences.add(out_labels[loc]);
        }
      }
    }
  }

  return relabel<OUT>(out_labels, voxels, next_label, equivalences, N, start_label);
}

template <typename LABEL, typename OUT = uint64_t>
OUT* connected_components(
  LABEL* in_labels, 
  const int64_t sx, const int64_t sy, const int64_t sz,
  std::vector<uint64_t> &num_components_per_slice,
  size_t &N = _dummy_N
) {

  const int64_t sxy = sx * sy;
  const int64_t voxels = sxy * sz;

  size_t max_labels = voxels;
  OUT* out_labels = new OUT[voxels]();
  N = 0;

  for (int64_t z = 0; z < sz; z++) {
    size_t tmp_N = 0;
    connected_components2d_4<LABEL, OUT>(
      (in_labels + sxy * z), sx, sy, 1, 
      max_labels, (out_labels + sxy * z),
      tmp_N, N
    );
    num_components_per_slice[z] = tmp_N;
    N += tmp_N;
  }

  return out_labels;
}


};

#endif
