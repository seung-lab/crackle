from typing import List

import numpy as np

from .headers import CrackleHeader, CrackFormat, LabelFormat
from . import crackcode
from .lib import width2dtype

def header(binary:bytes) -> dict:
  return CrackleHeader.frombytes(binary)

def labels(binary:bytes) -> np.ndarray:
  head = header(binary)

  if head.label_format == LabelFormat.FLAT:
    return np.unique(
      decode_flat_labels(binary, head.stored_dtype, head.dtype)
    )
  else:
    all_pins = raw_pins(binary)
    bgcolor = background_color(binary)
    return np.unique([ bgcolor ] + [ p['label'] for p in all_pins ])  

def raw_labels(binary:bytes) -> bytes:
  header = CrackleHeader.frombytes(binary)
  hb = header.HEADER_BYTES
  return binary[hb:hb+header.num_label_bytes]

def remap(binary:bytes, mapping:dict, preserve_missing_labels:bool = False):
  binary = bytearray(binary)
  header = CrackleHeader.frombytes(binary)
  hb = header.HEADER_BYTES
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

def components(binary:bytes):
  header = CrackleHeader.frombytes(binary)

  hl = len(header.tobytes())
  ll = header.num_label_bytes
  il = header.sz * header.z_index_width()
  cl = len(binary) - hl -ll - il

  return {
    'header': binary[:hl],
    'labels': binary[hl:hl+ll],
    'z_index': binary[hl+ll:hl+ll+il],
    'crack_codes': binary[-cl:],
  }

def component_lengths(binary:bytes):
  return { k:len(v) for k,v in components(binary).items() }

def background_color(binary:bytes) -> int:
  header = CrackleHeader.frombytes(binary)
  hb = CrackleHeader.HEADER_BYTES
  dtype = width2dtype[header.stored_data_width]
  bgcolor = np.frombuffer(binary[hb:hb+header.stored_data_width], dtype=dtype)
  return int(bgcolor[0])

def raw_pins(binary:bytes) -> np.ndarray:
  header = CrackleHeader.frombytes(binary)
  hb = CrackleHeader.HEADER_BYTES
  # background color followed by pins
  # but skip bgcolor
  pinset = binary[
    hb+header.stored_data_width:hb+header.num_label_bytes
  ]

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

  return [ 
    binary[code_offsets[i]:code_offsets[i+1]] 
    for i in range(header.sz) 
  ]

def crack_codes_to_cc_labels(
  crack_codes, 
  sx:int, sy:int, sz:int,
  permissible:bool
):
  cc_labels = np.zeros((sx,sy,sz), dtype=np.uint32)

  Ntotal = 0
  for z, code in enumerate(crack_codes):
    code = crackcode.unpack_binary(code, sx, sy)
    cc_slice, N = crackcode.decode_crack_code(
      code, sx, sy, permissible
    )
    cc_slice += Ntotal
    Ntotal += N
    cc_labels[:,:,z] = cc_slice

  return cc_labels, Ntotal

def decode_flat_labels(binary:bytes, stored_dtype, dtype):
  labels_binary = raw_labels(binary)
  return np.frombuffer(labels_binary, dtype=stored_dtype)\
    .astype(dtype, copy=False)

def decode_fixed_width_pins(
  binary:bytes, 
  cc_labels:np.ndarray, N:int, 
  label_dtype
):
  header = CrackleHeader.frombytes(binary)
  sx, sy, sz = header.sx, header.sy, header.sz

  bgcolor = background_color(binary)
  all_pins = raw_pins(binary)
  label_map = np.full((N,), fill_value=bgcolor, dtype=label_dtype)

  for pin in all_pins:
    for depth in range(pin['depth']+1):
      z = pin['idx'] // (sx*sy)
      y = (pin['idx'] - (z * sx*sy)) // sx
      x = pin['idx'] - (sx * y) - (sx * sy * z)
      ccid = cc_labels[x,y,z+depth]
      label_map[ccid] = pin['label']

  return label_map

def decompress(binary:bytes) -> np.ndarray:
  header = CrackleHeader.frombytes(binary)
  crack_codes = get_crack_codes(binary)

  sx,sy,sz = header.sx, header.sy, header.sz

  cc_labels, N = crack_codes_to_cc_labels(
    crack_codes, sx, sy, sz, 
    permissible=(header.crack_format == CrackFormat.PERMISSIBLE),
  )

  stored_label_dtype = width2dtype[header.stored_data_width]
  label_dtype = width2dtype[header.data_width]

  if header.label_format == LabelFormat.FLAT:
    label_map = decode_flat_labels(
      binary, stored_label_dtype, label_dtype
    )
  elif header.label_format == LabelFormat.PINS_FIXED_WIDTH:
    label_map = decode_fixed_width_pins(binary, cc_labels, N, label_dtype)
  else:
    raise ValueError(f"should never happen. labels fmt: {header.label_format}")

  out = np.zeros((sx,sy,sz), dtype=label_dtype)
  for z in range(sz):
    for y in range(sy):
      for x in range(sx):
        out[x,y,z] = label_map[cc_labels[x,y,z]]

  return out
