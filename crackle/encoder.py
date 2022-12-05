import numpy as np

import fastremap

from . import crackcode
from . import pins
from .header import CrackleHeader, LabelFormat, LabelSort, CrackFormat
from .lib import compute_byte_width, width2dtype
from .ccl import connected_components

def encode_flat_labels(labels, stored_data_dtype):
  sx,sy,sz = labels.shape

  cc_labels = np.zeros((sx,sy,sz), dtype=np.uint32)
  N_total = 0
  for z in range(sz):
    cc_slice, N = connected_components(labels[:,:,z])
    cc_slice += N_total
    cc_labels[:,:,z] = cc_slice
    N_total += N

  mapping = fastremap.component_map(cc_labels, labels)
  array = np.zeros((N_total,), dtype=stored_data_dtype)

  for ccid, label in mapping.items():
    array[ccid] = label

  return array.tobytes()

# parts of the file:
# HEADER, LABELS, ZINDEX, BOUNDARIES
def compress(labels:np.ndarray) -> bytes:
  if labels.ndim < 3:
    labels = labels[..., np.newaxis]

  sx,sy,sz = labels.shape

  num_pairs = fastremap.pixel_pairs(labels)
  crack_format = CrackFormat.IMPERMISSIBLE
  label_format = LabelFormat.PINS_FIXED_WIDTH
  if num_pairs / labels.size < 0.5:
    crack_format = CrackFormat.PERMISSIBLE
    label_format = LabelFormat.FLAT

  stored_data_width = compute_byte_width(np.max(labels))

  header = CrackleHeader(
    label_format=label_format,
    label_sort=LabelSort.INDEX,
    crack_format=crack_format,
    data_width=np.dtype(labels.dtype).itemsize,
    stored_data_width=stored_data_width,
    sx=labels.shape[0], sy=labels.shape[1], sz=labels.shape[2],
    num_label_bytes=0,
  )
  crack_codes = crackcode.encode_boundaries(
    labels, 
    permissible=(crack_format == CrackFormat.PERMISSIBLE),
  )
  z_index = [ len(code) for code in crack_codes ]

  if label_format == LabelFormat.PINS_FIXED_WIDTH:
    all_pins = pins.compute(labels)

    labels_binary = pins.fixed_width_binary(
      all_pins, 
      sx, sy, sz,
      stored_data_width=stored_data_width,
      index_width=header.index_width(),
      z_width=header.depth_width(),
      sort=header.label_sort,
    )
  elif label_format == LabelFormat.FLAT:
    labels_binary = encode_flat_labels(labels, width2dtype[stored_data_width])

  header.num_label_bytes = len(labels_binary)

  z_index_width = header.z_index_width()
  z_index = b''.join([ 
    sz.to_bytes(z_index_width, 'little') for sz in z_index 
  ])
  crack_codes = b''.join(crack_codes)

  return b''.join([
    header.tobytes(),
    labels_binary,
    z_index,
    crack_codes
  ])







