from typing import List

import numpy as np

from .header import CrackleHeader
from . import crackcode

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

def get_crack_codes(binary:bytes) -> List[bytes]:
  header = CrackleHeader.frombytes(binary)
  hb = CrackleHeader.HEADER_BYTES  
  offset = hb + header.num_label_bytes

  zindex_bytes = header.z_index_width() * header.sz

  z_index = np.frombuffer(
    binary[offset:offset+zindex_bytes], 
    dtype=width2dtype[header.z_index_width()]
  )

  code_offsets = np.concatenate((np.array([0]), z_index))
  code_offsets = np.add.accumulate(code_offsets)
  code_offsets += offset + zindex_bytes
  print(code_offsets)
  return [ 
    binary[code_offsets[i]:code_offsets[i+1]] 
    for i in range(header.sz) 
  ]

def crack_codes_to_cc_labels(
  crack_codes, 
  sx:int, sy:int, sz:int
):
  cc_labels = np.zeros((sx,sy,sz), dtype=np.uint32)

  Ntotal = 0
  for z, code in enumerate(crack_codes):
    code = crackcode.unpack_crack_binary(code)
    cc_slice, N = crackcode.decode_crack_code(code, sx, sy)
    cc_slice += Ntotal
    Ntotal += N
    cc_labels[:,:,z] = cc_slice

  return cc_labels, Ntotal

def decompress(binary: bytes) -> np.ndarray:
  header = CrackleHeader.frombytes(binary)
  hb = CrackleHeader.HEADER_BYTES
  all_pins = raw_pins(binary)
  crack_codes = get_crack_codes(binary)

  sx,sy,sz = header.sx, header.sy, header.sz

  cc_labels, N = crack_codes_to_cc_labels(
    crack_codes, sx, sy, sz
  )

  label_dtype = width2dtype[header.data_width]
  label_map = np.array((N,), dtype=label_dtype)

  for pin in all_pins:
    for depth in range(pin['depth']+1):
      z = pin['idx'] // (sx*sy)
      y = (pin['idx'] - (z * sx*sy)) // sx
      x = pin['idx'] - (sx * y) - (sx * sy * z)
      ccid = cc_labels[x,y,z+depth]
      label_map[ccid] = pin['label']

  out = np.zeros((sx,sy,sz), dtype=label_dtype)
  for z in range(sz):
    for y in range(sy):
      for x in range(sx):
        out[x,y,z] = label_map[cc_labels[x,y,z]]

  return out
