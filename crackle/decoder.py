import numpy as np

from .header import CrackleHeader

width2dtype = {
  1: np.uint8,
  2: np.uint16,
  4: np.uint32,
  8: np.uint64,
}

def labels(binary:bytes) -> np.ndarray:
  header = CrackleHeader.frombytes(binary)
  hb = CrackleHeader.HEADER_BYTES
  pinset = binary[hb:hb+header.num_label_bytes]
  dtype = np.dtype([
    ('label', width2dtype[header.stored_data_width]), 
    ('idx', width2dtype[header.index_width()]), 
    ('depth', width2dtype[header.depth_width()])
  ])
  all_pins = np.frombuffer(pinset, dtype=dtype)
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

def decompress(binary: bytes) -> np.ndarray:
  header = CrackleHeader.frombytes(binary)
  hb = CrackleHeader.HEADER_BYTES

  pinset = binary[hb:hb+header.num_label_bytes]
  dtype = np.dtype([
    ('label', width2dtype[header.stored_data_width]), 
    ('idx', width2dtype[header.index_width()]), 
    ('depth', width2dtype[header.depth_width()])
  ])
  all_pins = np.frombuffer(pinset, dtype=dtype)

  







