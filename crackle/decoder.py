from typing import List

import numpy as np

from .headers import CrackleHeader, CrackFormat, LabelFormat
from . import crackcode
from .lib import width2dtype, compute_byte_width, compute_dtype

def header(binary:bytes) -> dict:
  return CrackleHeader.frombytes(binary)

def labels(binary:bytes) -> np.ndarray:
  head = header(binary)

  if head.label_format == LabelFormat.FLAT:
    return np.unique(
      decode_flat_labels(binary, head.stored_dtype, head.dtype)
    )
  else:
    hb = CrackleHeader.HEADER_BYTES
    # bgcolor, num labels (u64), N labels, pins
    offset = hb + head.stored_data_width
    num_labels = int.from_bytes(binary[offset:offset+8], 'little')
    offset += 8
    labels = np.frombuffer(
      binary[offset:offset+num_labels*head.stored_data_width],
      dtype=head.stored_dtype
    )
    bgcolor = background_color(binary)
    labels = np.concatenate(([ bgcolor ], labels))
    labels.sort()
    return labels

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

def decode_pins(binary:bytes) -> np.ndarray:
  header = CrackleHeader.frombytes(binary)
  hb = CrackleHeader.HEADER_BYTES

  # bgcolor, num labels (u64), N labels, pins
  offset = hb + header.stored_data_width
  num_labels = int.from_bytes(binary[offset:offset+8], 'little')
  offset += 8
  labels = np.frombuffer(
    binary[offset:offset+num_labels*header.stored_data_width],
    dtype=header.stored_dtype
  )
  offset += num_labels*header.stored_data_width
  renum_width = compute_byte_width(len(labels))

  pinset = binary[
    offset:hb+header.num_label_bytes
  ]

  dtype = np.dtype([
    ('label', width2dtype[renum_width]),
    ('idx', width2dtype[header.index_width()]), 
    ('depth', width2dtype[header.depth_width()])
  ])
  pinset = np.frombuffer(pinset, dtype=dtype).copy()

  dtype_big = np.dtype([
    ('label', header.stored_dtype),
    ('idx', width2dtype[header.index_width()]), 
    ('depth', width2dtype[header.depth_width()])
  ])
  embiggened = np.empty((len(pinset),), dtype=dtype_big)
  embiggened[:] = pinset
  del pinset

  for pin in embiggened:
    pin['label'] = labels[pin['label']]
  return embiggened

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
  num_labels = int.from_bytes(labels_binary[:8], 'little')
  offset = 8

  uniq_bytes = num_labels * np.dtype(stored_dtype).itemsize
  uniq = np.frombuffer(labels_binary[8:8+uniq_bytes], dtype=stored_dtype)
  uniq = uniq.astype(dtype, copy=False)

  cc_label_dtype = compute_dtype(num_labels)
  cc_map = np.frombuffer(labels_binary[8+uniq_bytes:], dtype=cc_label_dtype)
  return uniq[cc_map]

def decode_fixed_width_pins(
  binary:bytes, 
  cc_labels:np.ndarray, N:int, 
  label_dtype
):
  header = CrackleHeader.frombytes(binary)
  sx, sy, sz = header.sx, header.sy, header.sz

  bgcolor = background_color(binary)
  all_pins = decode_pins(binary)
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
  elif header.label_format == LabelFormat.PINS_VARIABLE_WIDTH:
    raise NotImplemented("")
  else:
    raise ValueError(f"should never happen. labels fmt: {header.label_format}")

  out = np.zeros((sx,sy,sz), dtype=label_dtype)
  for z in range(sz):
    for y in range(sy):
      for x in range(sx):
        out[x,y,z] = label_map[cc_labels[x,y,z]]

  return out
