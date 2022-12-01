import numpy as np

from . import crackcode
from . import pins

class CrackleHeader:
  magic = b'crkl'
  format_version = 0

  def __init__(
    self, 
    data_width:int, stored_data_width:int,
    sx:int, sy:int, sz:int,
    num_label_bytes:int
  ):
    self.data_width = int(data_width)
    self.stored_data_width = int(stored_data_width)
    self.sx = int(sx)
    self.sy = int(sy)
    self.sz = int(sz)
    self.num_label_bytes = num_label_bytes
    # should we have a field that is y/n pins?

  @classmethod
  def frombytes(self, buffer:bytes):
    pass

  def tobytes(self) -> bytes:
    return b''.join([
      self.magic,
      self.format_version.to_bytes(1, 'little'),
      self.data_width.to_bytes(1, 'little'),
      self.stored_data_width.to_bytes(1, 'little'),
      self.sx.to_bytes(2, 'little'),
      self.sy.to_bytes(2, 'little'),
      self.sz.to_bytes(2, 'little'),
      self.num_label_bytes.to_bytes(8, 'little'),
    ])

def compute_byte_width(x) -> int:
  byte_width = 8
  if x < np.iinfo(np.uint8).max:
    byte_width = 1
  elif x < np.iinfo(np.uint16).max:
    byte_width = 2
  elif x < np.iinfo(np.uint32).max:
    byte_width = 4

  return byte_width

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
  z_levels = crackcode.encode_boundaries(labels)
  z_index = [ len(level) for level in z_levels ]

  all_pins = pins.compute(labels)
  labels_binary = pins.fixed_width_binary(
    all_pins, 
    sx, sy, sz,
    stored_data_width=stored_data_width,
    index_width=compute_byte_width(sx * sy * sz), 
    z_width=compute_byte_width(sx * sy),
  )

  header.num_label_bytes = len(labels_binary)

  z_index_width = compute_byte_width(sx * sy * 2)
  z_index = b''.join([ 
    sz.to_bytes(z_index_width, 'little') for sz in z_index 
  ])
  z_levels = b''.join(z_levels)

  return b''.join([
    header.tobytes(),
    labels_binary,
    z_index,
    z_levels
  ])







