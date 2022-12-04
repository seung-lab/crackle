from collections import defaultdict

import compresso
import numpy as np
import cc3d
from tqdm import tqdm

def extract_columns(labels:np.ndarray):
  sx,sy,sz = labels.shape

  cc_labels = np.zeros((sx,sy,sz), dtype=np.uint32)
  for z in range(sz):
    cc_labels[:,:,z] = cc3d.connected_components(labels[:,:,z], connectivity=4)

  pinsets = defaultdict(list)

  for x in range(sx):
    for y in range(sy):
      label = labels[x,y,0]
      label_set = set([ cc_labels[x,y,0] ])
      z_start = 0
      for z in range(1, sz):
        cur = labels[x,y,z]
        if label != cur:
          pinsets[label].append(
            ((x,y,z_start), (x,y,z), label_set)
          )
          label = cur
          z_start = z
          label_set = set()
        label_set |= set([ cc_labels[x,y,z] ])
      pinsets[label].append(
        ((x,y,z_start), (x,y,z), label_set)
      )

  return pinsets

def find_optimal_pins(pinset):
  sets = [ x[2] for x in pinset ]
  isets = [ [i,x] for i,x in enumerate(sets) ]

  universe = set.union(*sets)

  final_pins = []
  while len(universe):
    sizes = [ len(x[1]) for x in isets ]
    idx = sizes.index(max(sizes))
    i, cur = isets.pop(idx)
    universe -= cur
    for j, otherset in isets:
      otherset -= cur
    final_pins.append(pinset[i][:2])

  return final_pins

def compute(labels):
  pinsets = extract_columns(labels)

  all_pins = {}
  for label in pinsets:
    pin = find_optimal_pins(pinsets[label])
    all_pins[label] = pin

  return all_pins

def fixed_width_binary(
  all_pins, 
  sx:int, sy:int, sz:int,
  stored_data_width:int,
  index_width:int, z_width:int
) -> bytes:

  def toidx(tup):
    return int(tup[0] + sx * tup[1] + sx * sy * tup[2])

  linear = []
  for label, pins in all_pins.items():
    for pin in pins:
      linear.append(
        [ int(label), toidx(pin[0]), int(pin[1][2] - pin[0][2]) ]
      )

  linear = sorted(linear, key=lambda x: x[1])
  bytestream = []
  for pin in linear:
      bytestream.append(pin[0].to_bytes(stored_data_width, 'little'))
      bytestream.append(pin[1].to_bytes(index_width, 'little'))
      bytestream.append(pin[2].to_bytes(z_width, 'little'))
  del linear

  return b''.join(bytestream)










