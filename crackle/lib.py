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