import numpy as np
import google_crc32c

def crc32c(buffer) -> int:
  return int.from_bytes(
    google_crc32c.Checksum(buffer).digest(),
    'big'
  )

def compute_byte_width(x) -> int:
  byte_width = 8
  if x <= np.iinfo(np.uint8).max:
    byte_width = 1
  elif x <= np.iinfo(np.uint16).max:
    byte_width = 2
  elif x <= np.iinfo(np.uint32).max:
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

def eytzinger_sort_helper(inpt, output, i = 0, k = 1):
  """
  Takes an ascendingly sorted input and 
  an equal sized output buffer into which to 
  rewrite the input in eytzinger order.

  Modified from:
  https://algorithmica.org/en/eytzinger
  """
  if k <= len(inpt):
    i = eytzinger_sort_helper(inpt, output, i, 2 * k)
    output[k - 1] = inpt[i]
    i += 1
    i = eytzinger_sort_helper(inpt, output,i, 2 * k + 1)
  return i

def eytzinger_sort(labels):
  labels = sorted(labels)
  out = [0] * len(labels)
  eytzinger_sort_helper(labels, out)
  return out




