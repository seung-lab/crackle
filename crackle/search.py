from collections import defaultdict

import numpy as np

from .codec import header, decode_flat_labels

def compute_label_search_index(binary:bytes) -> Dict[List[int]]:
  """
  Scans the labels and computes a list of run pairs 
  [start, end) as half-open interval that describes all 
  grids needed to construct a label.

  This dict can be stored in a mapbuffer as a persistent index
  for decoding particular labels or finding their rough bounding
  box.
  """
  head = header(binary)
  elems = decode_flat_labels(head, binary)

  uniq = elems["unique"]
  cpg = elems["components_per_grid"]
  cpg = np.concatenate([ [0], cpg ])
  cpg = np.cumsum(cpg)
  keys = elems["cc_map"]

  index = defaultdict(list)

  for gi in range(head.num_grids()):
    uniq_keys = fastremap.unique(keys[cpg[gi]:cpg[gi+1]])
    for k in uniq_keys:
      seq = index[uniq[k]]
      if len(seq) and seq[-1] == gi:
        seq[-1] = gi+1
      else:
        seq.append(gi)
        seq.append(gi+1)
      
  return index