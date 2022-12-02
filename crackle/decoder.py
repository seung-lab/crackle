from typing import List

import numpy as np

from .header import CrackleHeader

width2dtype = {
  1: np.uint8,
  2: np.uint16,
  4: np.uint32,
  8: np.uint64,
}

def labels(binary:bytes) -> np.ndarray:
  all_pins = raw_pins(binary)
  return np.unique([ p['label'] for p in all_pins ])  

def remap(binary:bytes, mapping:dict, preserve_missing_labels:bool = False):
  binary = bytearray(binary)
  header = CrackleHeader.frombytes(binary)
  pinset = binary[hb:hb+header.num_label_bytes]
  dtype = np.dtype([
    ('label', width2dtype[header.stored_data_width]), 
    ('idx', width2dtype[header.index_width()]), 
    ('depth', width2dtype[header.depth_width()])
  ])
  all_pins = np.frombuffer(pinset, dtype=dtype)

  for pin in all_pins:
    try:
      ids[i] = mapping[ids[i]]
    except KeyError:
      if not preserve_missing_labels:
        raise

  return bytes(binary)

def nbytes(binary:bytes) -> np.ndarray:
  header = CrackleHeader.frombytes(binary)
  return header.data_width * header.sx * header.sy * header.sz

def raw_pins(binary:bytes) -> np.ndarray:
  header = CrackleHeader.frombytes(binary)
  hb = CrackleHeader.HEADER_BYTES
  pinset = binary[hb:hb+header.num_label_bytes]
  dtype = np.dtype([
    ('label', width2dtype[header.stored_data_width]), 
    ('idx', width2dtype[header.index_width()]), 
    ('depth', width2dtype[header.depth_width()])
  ])
  return np.frombuffer(pinset, dtype=dtype)

def get_crack(binary:bytes) -> List[bytes]:
  header = CrackleHeader.frombytes(binary)
  hb = CrackleHeader.HEADER_BYTES  
  offset = hb + header.num_label_bytes

  zindex_bytes = header.z_index_width() * header.sz
  z_index = np.frombuffer(
    binary[offset:offset+zindex_bytes], 
    dtype=width2dtype[header.z_index_width()]
  )
  # this calculation isn't right yet
  z_offsets = np.concatenate((np.array([0]), z_index))
  z_offsets = np.add.accumulate(z_index) + offset

  return [ binary[z_offsets[i]:z_offsets[i+1]] for i in range(sz) ]


def decompress(binary: bytes) -> np.ndarray:
  header = CrackleHeader.frombytes(binary)
  hb = CrackleHeader.HEADER_BYTES
  all_pins = raw_pins(binary)
  crack_codes = get_crack(binary)

  sx,sy,sz = header.sx, header.sy, header.sz

  # create CCL volume from crack codes
  # N = num CCL labels in volume
  cc_labels = np.zeros((header.sx, header.sy, header.sz)) # placeholder
  N = 100

  label_map = np.array((N,), dtype=width2dtype[header.data_width])

  for pin in all_pins:
    for depth in range(pin['depth']):
      z = pin['idx'] // (sx*sy)
      y = (pin['idx'] - (z * sx*sy)) // sx
      x = pin['idx'] - sx * (y - sy * z)
      ccid = cc_labels[x,y,z+depth]
      label_map[ccid] = pin['label']

  return label_map[cc_labels.flatten()].reshape((sx,sy,sz))




  







