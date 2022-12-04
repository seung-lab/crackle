import numpy as np

from . import crackcode
from . import pins
from .header import CrackleHeader, compute_byte_width

# parts of the file:
# HEADER, LABELS, ZINDEX, BOUNDARIES
def compress(labels:np.ndarray) -> bytes:
  if labels.ndim < 3:
    labels = labels[..., np.newaxis]

  sx,sy,sz = labels.shape

  stored_data_width = compute_byte_width(np.max(labels))

  header = CrackleHeader(
    data_width=np.dtype(labels.dtype).itemsize,
    stored_data_width=stored_data_width,
    sx=labels.shape[0], sy=labels.shape[1], sz=labels.shape[2],
    num_label_bytes=0,
  )
  crack_codes = crackcode.encode_boundaries(labels)
  z_index = [ len(code) for code in crack_codes ]

  all_pins = pins.compute(labels)
  labels_binary = pins.fixed_width_binary(
    all_pins, 
    sx, sy, sz,
    stored_data_width=stored_data_width,
    index_width=header.index_width(),
    z_width=header.depth_width(),
  )

  header.num_label_bytes = len(labels_binary)

  z_index_width = header.z_index_width()
  z_index = b''.join([ 
    sz.to_bytes(z_index_width, 'little') for sz in z_index 
  ])
  crack_codes = b''.join(crack_codes)

  print("header", header.tobytes())
  print("pins", labels_binary)
  print("z index", z_index)
  print("crack codes", crack_codes)

  return b''.join([
    header.tobytes(),
    labels_binary,
    z_index,
    crack_codes
  ])







