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

#ifndef __CC3D_SPECIAL_2_4_HPP__
#define __CC3D_SPECIAL_2_4_HPP__

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <stdexcept>

namespace crackle {
namespace cc3d {

static uint64_t _dummy_N;

template <typename T>
class DisjointSet {
public:
  T *ids;
  uint64_t length;

  DisjointSet (uint64_t len) {
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
      throw std::runtime_error("maximum length exception");
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
    uint64_t &N = _dummy_N, uint64_t start_label = 0
  ) {

  OUT label;
  std::unique_ptr<OUT[]> renumber(new OUT[num_labels + 1]());
  OUT next_label = start_label + 1;

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
  N = next_label - start_label - 1;
  for (int64_t loc = 0; loc < voxels; loc++) {
    out_labels[loc] = renumber[out_labels[loc]] - 1; // first label is 0 not 1
  }

  return out_labels;
}

uint64_t estimate_provisional_label_count_vcg(
  const std::vector<uint8_t>& vcg, const int64_t sx
) {
  uint64_t count = 0; // number of transitions between labels
  int64_t voxels = static_cast<int64_t>(vcg.size());

  for (int64_t loc = 0; loc < voxels; loc += sx) {
    count += 1;
    for (int64_t x = 1; x < sx; x++) {
      count += ((vcg[loc+x] & 0b0010) == 0);
    }
  }

  return count;
}

template <typename OUT>
OUT* color_connectivity_graph(
  const std::vector<uint8_t> &vcg, // voxel connectivity graph
  const int64_t sx, const int64_t sy, const int64_t sz,
  OUT* out_labels = NULL,
  uint64_t &N = _dummy_N
) {

  const int64_t sxy = sx * sy;
  const int64_t voxels = sx * sy * sz;

  uint64_t max_labels = estimate_provisional_label_count_vcg(vcg, sx) + 1;
  max_labels = std::min(max_labels, static_cast<uint64_t>(std::numeric_limits<OUT>::max()));

  if (out_labels == NULL) {
    out_labels = new OUT[voxels]();
  }

  DisjointSet<OUT> equivalences(max_labels);

  OUT new_label = 0;
  for (int64_t z = 0; z < sz; z++) {
    new_label++;
    equivalences.add(new_label);

    for (int64_t x = 0; x < sx; x++) {
      if (x > 0 && (vcg[x + sxy * z] & 0b0010) == 0) {
        new_label++;
        equivalences.add(new_label);
      }
      out_labels[x + sxy * z] = new_label;
    }

    const int64_t B = -1;
    const int64_t C = -sx;

    for (int64_t y = 1; y < sy; y++) {
      for (int64_t x = 0; x < sx; x++) {
        int64_t loc = x + sx * y + sxy * z;

        if (x > 0 && (vcg[loc] & 0b0010)) {
          out_labels[loc] = out_labels[loc+B];
          if (y > 0 && (vcg[loc + C] & 0b0010) == 0 && (vcg[loc] & 0b1000)) {
            equivalences.unify(out_labels[loc], out_labels[loc+C]);
          }
        }
        else if (y > 0 && vcg[loc] & 0b1000) {
          out_labels[loc] = out_labels[loc+C];
        }
        else {
          new_label++;
          out_labels[loc] = new_label;
          equivalences.add(new_label);
        }
      }
    }
  }

  relabel<OUT>(out_labels, voxels, new_label, equivalences, N);
  return out_labels;
}

template <typename LABEL>
uint64_t estimate_provisional_label_count(
  const LABEL* in_labels, const int64_t sx, const int64_t voxels
) {
  uint64_t count = 0; // number of transitions between labels

  for (int64_t loc = 0; loc < voxels; loc += sx) {
    count += 1;
    for (int64_t x = 1; x < sx; x++) {
      count += (in_labels[loc + x] != in_labels[loc + x - 1]);
    }
  }

  return count;
}

template <typename LABEL, typename OUT>
OUT* connected_components2d_4(
    const LABEL* in_labels, 
    const int64_t sx, const int64_t sy, const int64_t sz,
    OUT *out_labels = NULL, 
    const uint64_t start_label = 0, uint64_t &N = _dummy_N
  ) {

  const int64_t sxy = sx * sy;
  const int64_t voxels = sx * sy * sz;

  uint64_t max_labels = estimate_provisional_label_count<LABEL>(in_labels, sx, voxels) + 1;
  max_labels = std::min(max_labels, static_cast<uint64_t>(std::numeric_limits<OUT>::max()));

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
  const LABEL* in_labels, 
  const int64_t sx, const int64_t sy, const int64_t sz,
  std::vector<uint64_t> &num_components_per_slice,
  OUT* out_labels = NULL,
  uint64_t &N = _dummy_N
) {

  const int64_t sxy = sx * sy;
  N = 0;

  if (out_labels == NULL) {
    out_labels = new OUT[sx * sy * sz]();
  }

  for (int64_t z = 0; z < sz; z++) {
    uint64_t tmp_N = 0;
    connected_components2d_4<LABEL, OUT>(
      (in_labels + sxy * z), 
      sx, sy, 1, 
      (out_labels + sxy * z),
       N, tmp_N
    );
    num_components_per_slice[z] = tmp_N;
    N += tmp_N;
  }

  return out_labels;
}


};
};

#endif
