from collections import defaultdict

import numpy as np
from tqdm import tqdm

import fastcrackle

from .lib import compute_byte_width, eytzinger_sort

def append_pin(pinsets, label, z_start, x, y, z, cc_set):
  # try to reduce candidate pins by filtering
  # out ones that are strictly dominated by their
  # neighboring pin
  if (
    len(pinsets[label])
    and pinsets[label][-1][0][0] == x-1
    and pinsets[label][-1][0][1] == y
  ):
    if (
      pinsets[label][-1][0][2] <= z_start
      and pinsets[label][-1][1][2] >= z
    ):
      pass
    elif (
      pinsets[label][-1][0][2] >= z_start
      and pinsets[label][-1][1][2] <= z
    ):
      pinsets[label][-1] = (
        (x,y,z_start), (x,y,z), set(cc_set)
      )
    else:
      pinsets[label].append(
        ((x,y,z_start), (x,y,z), set(cc_set))
      )
  else:
    pinsets[label].append(
      ((x,y,z_start), (x,y,z), set(cc_set))
    )

def extract_columns(labels:np.ndarray):
  sx,sy,sz = labels.shape

  cc_labels, components, N_total = fastcrackle.connected_components(labels)
  cc_labels = cc_labels.reshape(labels.shape, order='F')

  pinsets = defaultdict(list)
    
  for y in range(sy):
    for x in range(sx):
      label = labels[x,y,0]
      z_start = 0
      for z in range(1, sz):
        cur = labels[x,y,z]
        if label != cur:
          append_pin(pinsets, label, z_start, x, y, z-1, cc_labels[x,y,z_start:z])
          label = cur
          z_start = z

      append_pin(pinsets, label, z_start, x, y, z, cc_labels[x,y,z_start:])

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
    isets = [ x for x in isets if len(x[1]) > 0 ]
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
  """Format is [label][pin top index][pin depth] for each pin."""

  def toidx(tup):
    return int(tup[0] + sx * tup[1] + sx * sy * tup[2])

  # find bg color
  max_pins_label = 0
  max_pins = 0
  for label, pins in all_pins.items():
    if len(pins) > max_pins:
      max_pins_label = int(label)
      max_pins = len(pins)

  all_pins.pop(max_pins_label)

  linear = []
  for label, pins in all_pins.items():
    for pin in pins:
      linear.append(
        [ int(label), toidx(pin[0]), int(pin[1][2] - pin[0][2]) ]
      )

  linear = sorted(linear, key=lambda x: x[1])
  bgcolor = max_pins_label.to_bytes(stored_data_width, 'little')

  all_labels = eytzinger_sort(list(all_pins.keys()))

  bytestream = []
  bytestream.append(bgcolor) # bgcolor
  bytestream.append(len(all_labels).to_bytes(8, 'little'))
  for label in all_labels:
    bytestream.append(int(label).to_bytes(stored_data_width, 'little'))

  renumbering = { x: i for i, x in enumerate(all_labels) }
  renum_data_width = compute_byte_width(len(renumbering))

  for pin in linear:
    bytestream.append(
      renumbering[pin[0]].to_bytes(renum_data_width, 'little')
    )
    bytestream.append(pin[1].to_bytes(index_width, 'little'))
    bytestream.append(pin[2].to_bytes(z_width, 'little'))
  del linear

  return b''.join(bytestream)

def condensed_binary(
  all_pins, 
  sx:int, sy:int, sz:int,
  stored_data_width:int,
  index_width:int, z_width:int
):
  """Format is [label][num pins][pin1]...[pin_N] for each label."""

  def toidx(tup):
    return int(tup[0] + sx * tup[1] + sx * sy * tup[2])

  linear = []
  for label, pins in all_pins.items():
    seq = []
    seq.append(int(label))
    seq.append(len(pins))
    for pin in pins:
      seq.append(
        [ toidx(pin[0]), int(pin[1][2] - pin[0][2]) ]
      )
    linear.append(seq)

  bytestream = []
  for seq in linear:
    bytestream.append(seq[0].to_bytes(stored_data_width, 'little'))
    bytestream.append(seq[1].to_bytes(1, 'little'))
    for pin in seq[2:]:
      bytestream.append(pin[0].to_bytes(index_width, 'little'))
      bytestream.append(pin[1].to_bytes(z_width, 'little'))
  del linear

  return b''.join(bytestream)
