from typing import List, Optional
from collections import namedtuple

import numpy as np
import fastremap
import fastcrackle

from .headers import CrackleHeader, CrackFormat, LabelFormat, FormatError
from .lib import width2dtype, compute_byte_width, compute_dtype

def header(binary:bytes) -> CrackleHeader:
  """Decode the header from a Crackle bytestream."""
  return CrackleHeader.frombytes(binary)

def labels(binary:bytes) -> np.ndarray:
  """Extract the unique labels from a Crackle bytestream."""
  head = header(binary)
  if head.voxels() == 0:
    return np.zeros((0,), dtype=head.dtype)

  hb = CrackleHeader.HEADER_BYTES
  offset = hb + head.sz * 4

  if head.label_format == LabelFormat.FLAT:
    # num labels (u64), N labels
    num_labels = int.from_bytes(binary[offset:offset+8], 'little')
    offset += 8
    return np.frombuffer(
      binary[offset:offset+num_labels*head.stored_data_width],
      dtype=head.stored_dtype
    ).astype(head.dtype, copy=False)
  else:
    # bgcolor, num labels (u64), N labels, pins
    offset += head.stored_data_width
    num_labels = int.from_bytes(binary[offset:offset+8], 'little')
    offset += 8
    labels = np.frombuffer(
      binary[offset:offset+num_labels*head.stored_data_width],
      dtype=head.stored_dtype
    )
    bgcolor = background_color(binary)
    labels = np.concatenate(([ bgcolor ], labels))
    labels.sort()
    return labels.astype(head.dtype, copy=False)

def contains(binary:bytes, label:int) -> bool:
  """Rapidly check if a label exists in a Crackle bytestream."""
  head = header(binary)
  hb = CrackleHeader.HEADER_BYTES
  offset = hb + head.sz * 4

  # bgcolor, num labels (u64), N labels, pins
  if head.label_format == LabelFormat.PINS_FIXED_WIDTH:
    bgcolor = background_color(binary)
    if bgcolor == label:
      return True
    offset += head.stored_data_width

  num_labels = int.from_bytes(binary[offset:offset+8], 'little')
  offset += 8
  uniq = np.frombuffer(
    binary[offset:offset+num_labels*head.stored_data_width],
    dtype=head.stored_dtype
  )
  idx = np.searchsorted(uniq, label)
  if 0 < idx < uniq.size:
    return True
  elif idx == 0:
    return uniq[0] == label
  else:
    return False

def raw_labels(binary:bytes) -> bytes:
  header = CrackleHeader.frombytes(binary)
  offset = header.HEADER_BYTES + header.sz * 4
  return binary[offset:offset+header.num_label_bytes]

def remap(binary:bytes, mapping:dict, preserve_missing_labels:bool = False):
  """
  Remap the labels in a crackle bystream without decompressing.
  
  binary: A Crackle bytestream
  mapping: A dict mapping old labels to new labels.
  preserve_missing_labels: By default we presume all values
    will be remapped. If a value in the binary does not have
    a corresponding mapping key, it will raise a KeyError. If
    True, just leave it be and encode all the labels that are 
    in the mapping.
  """
  orig = binary
  binary = bytearray(binary)
  
  head = CrackleHeader.frombytes(binary)
  hb = CrackleHeader.HEADER_BYTES

  # flat: num_labels, N labels, remapped labels
  # pins: bgcolor, num labels (u64), N labels, pins

  offset = hb + 4 * head.sz
  if head.label_format == LabelFormat.PINS_FIXED_WIDTH:
    bgcolor = int.from_bytes(binary[offset:offset+head.stored_data_width], 'little')

    try:
      binary[offset:offset+head.stored_data_width] = \
        mapping[bgcolor].to_bytes(head.stored_data_width, 'little')
    except KeyError:
      if not preserve_missing_labels:
        raise

    offset += head.stored_data_width

  num_labels = int.from_bytes(binary[offset:offset+8], 'little')
  offset += 8
  uniq_bytes = num_labels * head.stored_data_width
  all_labels = np.frombuffer(
    binary[offset:offset+uniq_bytes],
    dtype=head.stored_dtype
  )
  fastremap.remap(
    all_labels, mapping, 
    preserve_missing_labels=preserve_missing_labels, 
    in_place=True
  )
  binary[offset:offset+uniq_bytes] = list(all_labels.view(np.uint8))
  return bytes(binary)

def refit(binary:bytes) -> bytes:
  """
  Change the rendered dtype to the smallest
  dtype needed to render the image without
  loss of precision.
  """
  head = header(binary)
  head.data_width = head.stored_data_width
  return b''.join([ 
    head.tobytes(), 
    binary[CrackleHeader.HEADER_BYTES:] 
  ])

def nbytes(binary:bytes) -> np.ndarray:
  """Compute the size in bytes of the decompressed array."""
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
    'z_index': binary[hl:hl+il],
    'labels': binary[hl+il:hl+ll+il],
    'crack_codes': binary[-cl:],
  }

def component_lengths(binary:bytes):
  return { k:len(v) for k,v in components(binary).items() }

def crack_codes(binary:bytes) -> np.ndarray:
  header = CrackleHeader.frombytes(binary)
  comps = components(binary)
  z_index = np.frombuffer(comps["z_index"], dtype=np.uint32)
  z_index = np.cumsum(z_index)
  z_index += (
    len(header.tobytes()) 
    + header.num_label_bytes 
    + header.sz * header.z_index_width()
  )
  z_index = np.concatenate((z_index, [ len(binary) ]))
  z_index = z_index.astype(np.uint64)
  
  codes = []
  for i in range(header.sz):
    codes.append(
      binary[z_index[i]:z_index[i+1]]
    )
  return codes

def boc(crack_codes:bytes) -> np.ndarray:
  """extract the beginning of chain region from the crack code"""
  N = int.from_bytes(crack_codes[:4], byteorder='little')
  return crack_codes[:N+4]

def background_color(binary:bytes) -> int:
  """For pin encodings only, extract the background color."""
  header = CrackleHeader.frombytes(binary)

  if header.label_format == LabelFormat.FLAT:
    raise FormatError("Background color can only be extracted from pin encoded streams.")

  offset = CrackleHeader.HEADER_BYTES + header.sz * 4
  dtype = width2dtype[header.stored_data_width]
  bgcolor = np.frombuffer(binary[offset:offset+header.stored_data_width], dtype=dtype)
  return int(bgcolor[0])

def decode_pins(binary:bytes) -> np.ndarray:
  header = CrackleHeader.frombytes(binary)

  if header.label_format == LabelFormat.PINS_FIXED_WIDTH:
    return decode_fixed_pins(binary)
  elif header.label_format == LabelFormat.PINS_VARIABLE_WIDTH:
    return decode_condensed_pins(binary)
  else:
    raise FormatError("Cannot decode pins from flat format.")

def decode_condensed_pins(binary:bytes) -> np.ndarray:
  header = CrackleHeader.frombytes(binary)
  hb = CrackleHeader.HEADER_BYTES

  if header.label_format != LabelFormat.PINS_VARIABLE_WIDTH:
    raise FormatError("This function can only extract pins from variable width streams.")

  # bgcolor, num labels (u64), N labels, pins
  offset = hb + header.sz * 4 + header.stored_data_width
  labels_end = hb + header.sz * 4 + header.num_label_bytes
  num_labels = int.from_bytes(binary[offset:offset+8], 'little')
  offset += 8
  uniq = np.frombuffer(
    binary[offset:offset+num_labels*header.stored_data_width],
    dtype=header.stored_dtype
  )
  offset += num_labels * header.stored_data_width  
  combined_width = binary[offset]

  num_pins_width = 2 ** (combined_width & 0b11)
  depth_width = 2 ** ((combined_width >> 2) & 0b11)
  index_width = header.index_width()

  offset += 1
  pinset = binary[offset:labels_end]
  idtype = np.dtype(width2dtype[header.index_width()])
  ddtype = np.dtype(width2dtype[depth_width])

  PinTuple = namedtuple('Pin', ['index', 'depth'])

  pins = {}

  offset = 0
  for label in range(num_labels):
    n_pins = int.from_bytes(pinset[offset:offset+num_pins_width], 'little')
    offset += num_pins_width
    index_arr = np.frombuffer(pinset[offset:offset+n_pins*idtype.itemsize], dtype=idtype)
    index_arr = index_arr.copy()
    for i in range(1, len(index_arr)):
      index_arr[i] += index_arr[i-1]
    offset += n_pins*idtype.itemsize
    depth_arr = np.frombuffer(pinset[offset:offset+n_pins*ddtype.itemsize], dtype=ddtype)
    offset += n_pins * ddtype.itemsize
    pins[uniq[label]] = [ PinTuple(i,d) for i,d in zip(index_arr, depth_arr) ]
    
  return pins

def decode_fixed_pins(binary:bytes) -> np.ndarray:
  """For pin encodings only, extract the pins."""
  header = CrackleHeader.frombytes(binary)
  hb = CrackleHeader.HEADER_BYTES

  if header.label_format != LabelFormat.PINS_FIXED_WIDTH:
    raise FormatError("This function can only extract pins from fixed width streams.")

  # bgcolor, num labels (u64), N labels, pins
  offset = hb + header.stored_data_width + header.sz * 4
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

def decode_flat_labels(binary:bytes, stored_dtype, dtype, sz:int):
  labels_binary = raw_labels(binary)
  num_labels = int.from_bytes(labels_binary[:8], 'little')
  offset = 8

  uniq_bytes = num_labels * np.dtype(stored_dtype).itemsize
  uniq = np.frombuffer(labels_binary[8:8+uniq_bytes], dtype=stored_dtype)
  uniq = uniq.astype(dtype, copy=False)

  cc_label_dtype = compute_dtype(num_labels)
  cc_map = np.frombuffer(labels_binary[8+uniq_bytes+4*sz:], dtype=cc_label_dtype)
  return uniq[cc_map]

def decompress(binary:bytes) -> np.ndarray:
  """Decompress a Crackle binary into a Numpy array."""
  return decompress_range(binary, None, None)

def decompress_range(binary:bytes, z_start:Optional[int], z_end:Optional[int]) -> np.ndarray:
  """
  Decompress a Crackle binary into a Numpy array.

  A partial z-range can be decoded by setting z_start and z_end.
  """
  header = CrackleHeader.frombytes(binary)
  sx, sy, sz = header.sx, header.sy, header.sz

  if z_start is None:
    z_start = 0
  if z_end is None:
    z_end = sz

  if (sx * sy * sz == 0):
    labels = np.zeros((0,), dtype=header.dtype)
  else:
    labels = fastcrackle.decompress(binary, z_start, z_end)

  order = 'F' if header.fortran_order else 'C'
  labels = labels.reshape((sx,sy,z_end - z_start), order=order)

  if header.signed:
    if header.data_width == 1:
      labels = labels.view(np.int8)
    elif header.data_width == 2:
      labels = labels.view(np.int16)
    elif header.data_width == 4:
      labels = labels.view(np.int32)
    else:
      labels = labels.view(np.int64)

  return labels

def compress(
  labels:np.ndarray, 
  allow_pins:bool = False,
  markov_model_order:int = 0
) -> bytes:
  """
  Compress the 3D labels array into a Crackle bytestream.

  [EXPERIMENTAL]
  BINARIES ENCODED USING THIS OPTION MAY BREAK IN FUTURE VERSIONS
  allow_pins: If True, when the number of voxel pairs in the 
    volume is > 50% of the number of voxels, use the pin encoding
    strategy. Pins use 3D information to encode the label map.
    However, they are still investigational and currently work 
    best on simpler images. Pin computation requires appoximately
    solving a set cover problem and can be slow on larger images.
    However, it can likely be improved.
  """
  if np.issubdtype(labels.dtype, np.signedinteger):
    raise TypeError("Signed integer data types are not currently supported.")

  f_order = labels.flags.f_contiguous
  labels = np.asfortranarray(labels)
  return fastcrackle.compress(labels, allow_pins, f_order, markov_model_order)

