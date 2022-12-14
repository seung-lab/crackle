import numpy as np

def compute_byte_width(x) -> int:
  byte_width = 8
  if x < np.iinfo(np.uint8).max:
    byte_width = 1
  elif x < np.iinfo(np.uint16).max:
    byte_width = 2
  elif x < np.iinfo(np.uint32).max:
    byte_width = 4

  return byte_width

width2dtype = {
  1: np.uint8,
  2: np.uint16,
  4: np.uint32,
  8: np.uint64,
}

def compute_dtype(x) -> np.dtype:
  return width2dtype[compute_byte_width(x)]

def pack_bits(attrs):
  """[ (3, 2), (12, 4) ] # (value, bits)"""
  encoded = 0
  offset = 0
  for value, bits in attrs:
    assert value < (2 ** bits)
    encoded = encoded | (int(value) << offset)
    offset += bits
  return encoded

def unpack_bits(encoded, bits_per_value):
  """[1,3,2,4... etc]"""
  unpacked = []
  offset = 0
  for bits in bits_per_value:
    unpacked.append(
      (encoded >> offset) & ((1 << bits) - 1)
    )
    offset += bits
  return unpacked





