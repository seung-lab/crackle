from typing import List, Optional, Tuple, Sequence, Union
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

def num_labels(binary:bytes) -> int:
  """Returns the number of unique labels."""
  head = header(binary)
  hb = CrackleHeader.HEADER_BYTES

  if head.voxels() == 0:
    return 0

  offset = hb + head.sz * 4
  N = 0
  if head.label_format != LabelFormat.FLAT:
    offset += head.stored_data_width
    N += 1 # bgcolor
  N += int.from_bytes(binary[offset:offset+8], 'little')
  return N

def contains(binary:bytes, label:int) -> bool:
  """Rapidly check if a label exists in a Crackle bytestream."""
  head = header(binary)

  if not head.is_sorted:
    return label in labels(binary)

  hb = CrackleHeader.HEADER_BYTES
  offset = hb + head.sz * 4

  # bgcolor, num labels (u64), N labels, pins
  if head.label_format == LabelFormat.PINS_VARIABLE_WIDTH:
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
  if 0 <= idx < uniq.size:
    return uniq[idx] == label
  else:
    return False

def raw_labels(binary:bytes) -> bytes:
  header = CrackleHeader.frombytes(binary)
  offset = header.HEADER_BYTES + header.sz * 4
  return binary[offset:offset+header.num_label_bytes]

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
  z_index = np.concatenate([ [0], z_index ])
  z_index = np.cumsum(z_index)
  z_index += (
    len(header.tobytes()) 
    + header.num_label_bytes 
    + header.sz * header.z_index_width()
  )
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
    return decode_condensed_pins(binary)[0]
  else:
    raise FormatError("Cannot decode pins from flat format.")

def decode_condensed_pins_components(binary:bytes) -> dict:
  components = {}
  head = CrackleHeader.frombytes(binary)
  hb = CrackleHeader.HEADER_BYTES

  if head.label_format != LabelFormat.PINS_VARIABLE_WIDTH:
    raise FormatError("This function can only extract pins from variable width streams.")

  # bgcolor, num labels (u64), N labels, pins
  labels_binary = raw_labels(binary)
  bgcolor = background_color(binary)
  offset = head.stored_data_width
  num_labels = int.from_bytes(labels_binary[offset:offset+8], 'little')
  offset += 8
  uniq = np.frombuffer(
    labels_binary[offset:offset+num_labels*head.stored_data_width],
    dtype=head.stored_dtype
  )
  offset += num_labels * head.stored_data_width  
  
  component_dtype = width2dtype[head.component_width()]
  component_bytes = head.num_grids() * head.component_width()
  components_per_grid = np.frombuffer(
    labels_binary[offset:offset+component_bytes], 
    dtype=component_dtype
  )
  offset += component_bytes

  combined_width = labels_binary[offset]
  offset += 1

  num_pins_width = 2 ** (combined_width & 0b11)
  depth_width = 2 ** ((combined_width >> 2) & 0b11)
  cc_labels_width = 2 ** ((combined_width >> 4) & 0b11)

  pinset = labels_binary[offset:]

  return {
    "bgcolor": bgcolor,
    "uniq": uniq,
    "components_per_grid": components_per_grid,
    "num_pins_width": num_pins_width,
    "depth_width": depth_width,
    "cc_labels_width": cc_labels_width,
    "pinset": pinset,
  }

def decode_condensed_pins(binary:bytes) -> np.ndarray:
  head = CrackleHeader.frombytes(binary)

  if head.label_format != LabelFormat.PINS_VARIABLE_WIDTH:
    raise FormatError("This function can only extract pins from variable width streams.")

  elems = decode_condensed_pins_components(binary)
  components_per_grid = elems["components_per_grid"]
  components_per_grid = np.cumsum(components_per_grid)

  num_pins_width = elems["num_pins_width"]
  depth_width = elems["depth_width"]
  cc_labels_width = elems["cc_labels_width"]
  uniq = elems["uniq"]

  pinset = elems["pinset"]

  idtype = np.dtype(width2dtype[head.index_width()])
  ddtype = np.dtype(width2dtype[depth_width])
  cdtype = np.dtype(width2dtype[cc_labels_width])

  PinTuple = namedtuple('Pin', ['index', 'depth'])

  pins = {}
  single_labels = {}

  offset = 0
  for label in range(uniq.size):
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

    num_cc_labels = int.from_bytes(pinset[offset:offset+num_pins_width], 'little')
    offset += num_pins_width
    cc_labels = np.frombuffer(pinset[offset:offset+num_cc_labels*cc_labels_width], dtype=cdtype)
    cc_labels = np.cumsum(cc_labels)
    offset += num_cc_labels * cc_labels_width

    single_labels[uniq[label]] = cc_labels

  return pins, single_labels

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
  head = header(binary)
  if head.label_format != LabelFormat.FLAT:
    raise FormatError("Must be flat labels format.")

  labels_binary = raw_labels(binary)
  num_labels = int.from_bytes(labels_binary[:8], 'little')
  offset = 8

  uniq_bytes = num_labels * np.dtype(stored_dtype).itemsize
  uniq = np.frombuffer(labels_binary[offset:offset+uniq_bytes], dtype=stored_dtype)
  uniq = uniq.astype(dtype, copy=False)

  offset += uniq_bytes
  component_dtype = width2dtype[head.component_width()]
  component_bytes = head.num_grids() * head.component_width()
  components_per_grid = np.frombuffer(
    labels_binary[offset:offset+component_bytes], 
    dtype=component_dtype
  )

  offset += component_bytes

  cc_label_dtype = compute_dtype(num_labels)
  cc_map = np.frombuffer(labels_binary[offset:], dtype=cc_label_dtype)
  
  return {
    "num_labels": num_labels,
    "unique": uniq,
    "components_per_grid": components_per_grid,
    "cc_map": cc_map,
  }

def z_range_for_label(binary:bytes, label:int) -> Tuple[int,int]:
  head = header(binary)
  if head.label_format == LabelFormat.FLAT:
    return z_range_for_label_flat(binary, label)
  elif head.label_format == LabelFormat.PINS_VARIABLE_WIDTH:
    return z_range_for_label_condensed_pins(binary, label)
  else:
    raise ValueError("Label format not supported.")

def z_range_for_label_flat(binary:bytes, label:int) -> Tuple[int,int]:
  head = header(binary)
  labels_binary = raw_labels(binary)
 
  num_labels = int.from_bytes(labels_binary[:8], 'little')
  offset = 8
  uniq = np.frombuffer(
    labels_binary[offset:offset+num_labels*head.stored_data_width],
    dtype=head.stored_dtype
  )
  idx = np.searchsorted(uniq, label)
  if idx < 0 or idx >= uniq.size or uniq[idx] != label:
    return (-1, -1)

  offset += num_labels * head.stored_data_width
  next_offset = offset + head.num_grids() * head.component_width()
  dtype = width2dtype[head.component_width()]

  components_per_grid = np.frombuffer(labels_binary[offset:next_offset], dtype=dtype)
  components_per_grid = np.cumsum(components_per_grid)

  offset = next_offset

  dtype = compute_dtype(num_labels)
  cc_labels = np.frombuffer(labels_binary[offset:], dtype=dtype)

  cc_idxs = np.where(cc_labels == idx)[0]

  if cc_idxs.size == 0:
    return (-1, -1)

  min_cc = cc_idxs[0]
  max_cc = cc_idxs[-1]

  z_start = 0
  z_end = head.sz - 1

  for z in range(head.sz):
    if components_per_grid[z] >= min_cc:
      z_start = z
      break

  for z in range(head.sz - 1, -1, -1):
    if components_per_grid[z] <= max_cc:
      z_end = z + 1
      break

  return (z_start, z_end+1)

def z_range_for_label_condensed_pins(binary:bytes, label:int) -> Tuple[int,int]:
  head = header(binary)
  labels_binary = raw_labels(binary)

  bgcolor = background_color(binary)
  if bgcolor == label:
    return (0, head.sz)
  
  offset = head.stored_data_width
  num_labels = int.from_bytes(labels_binary[offset:offset+8], 'little')
  offset += 8
  uniq = np.frombuffer(
    labels_binary[offset:offset+num_labels*head.stored_data_width],
    dtype=head.stored_dtype
  )
  idx = np.searchsorted(uniq, label)
  if idx < 0 or idx >= uniq.size or uniq[idx] != label:
    return (-1, -1)

  offset += num_labels*head.stored_data_width
  component_dtype = width2dtype[head.component_width()]
  component_bytes = head.num_grids() * head.component_width()
  components_per_grid = np.frombuffer(
    labels_binary[offset:offset+component_bytes], 
    dtype=component_dtype
  )
  components_per_grid = np.cumsum(components_per_grid)
  all_pins, all_single_labels = decode_condensed_pins(binary)
  label_pins = all_pins[label]
  single_labels = all_single_labels[label]

  z_start = head.sz - 1
  z_end = 0

  sxy = head.sx * head.sy
  min = __builtins__["min"]
  max = __builtins__["max"]

  for pin in label_pins:
    z = pin.index // sxy
    z_start = min(z_start, z)
    z_end = max(z_end, z+pin.depth + 1)

  if len(single_labels) == 0:
    return (z_start, z_end)

  for lbl in  [ single_labels[0], single_labels[-1] ]:
    z = np.searchsorted(components_per_grid, lbl) - 1
    z_start = min(z_start, z)
    z_end = max(z_end, z)

  z_start = max(z_start, 0)
  z_end = min(z_end + 2, head.sz)

  return (z_start, z_end)

def decompress_binary_image(binary:bytes, label:Optional[int]) -> np.ndarray:
  z_start, z_end = z_range_for_label(binary, label)
  header = CrackleHeader.frombytes(binary)
  order = "F" if header.fortran_order else "C"
  image = np.zeros([header.sx, header.sy, header.sz], dtype=bool, order=order)

  if z_start == -1 and z_end == -1:
    return image

  cutout = decompress_range(binary, z_start, z_end)
  image[:,:,z_start:z_end] = (cutout == label)
  return image

def decompress(binary:bytes, label:Optional[int] = None) -> np.ndarray:
  """
  Decompress a Crackle binary into a Numpy array. 
  If label is provided, decompress into  a binary (bool) image.
  """
  if label is None:
    return decompress_range(binary, None, None)
  return decompress_binary_image(binary, label)

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

def extract_keys(binary:bytes) -> np.ndarray:
  head = header(binary)
  if head.label_format != LabelFormat.FLAT:
    raise FormatError("Can't use this function except with FLAT labels.")

  N = num_labels(binary)
  raw = raw_labels(binary)
  idx_bytes = head.component_width() * head.sz
  offset = 8 + N * head.stored_data_width + idx_bytes
  key_width = compute_byte_width(N)
  return np.frombuffer(raw[offset:], dtype=f'u{key_width}')

def condense_unique(binary:bytes) -> bytes:
  """
  A remapped crackle array may have
  elements in its unique array that are
  duplicated or not sorted. This will
  condense the information in the array
  and set the "is_sorted" flat to True.

  Note that fully decompressing and recompressing
  may still yield benefits as the crack code
  will not be adjusted. If two adjacent 
  connected components were mapped to the
  same label, the crack code will remain 
  oversegmented.
  """
  head = header(binary)
  uniq = labels(binary)

  reduced_uniq = fastremap.unique(uniq)

  if np.all(uniq == reduced_uniq):
    return binary

  mapping = { u: i for i, u in enumerate(reduced_uniq) }

  N = num_labels(binary)
  key_width = compute_byte_width(N)

  keys = extract_keys(binary)
  keys = np.array([ mapping[u] for u in uniq[keys] ], dtype=f'u{key_width}')

  raw = raw_labels(binary)
  idx_bytes = head.component_width() * head.sz
  offset = 8 + N * head.stored_data_width + idx_bytes
  labels_idx = raw[offset:offset + idx_bytes]

  labels_binary = b''.join([
    len(reduced_uniq).to_bytes(8, 'little'),
    reduced_uniq.astype(head.stored_dtype).tobytes(),
    labels_idx,
    keys.tobytes(),
  ])

  comps = components(binary)
  head.num_label_bytes = len(labels_binary)
  head.is_sorted = True

  return b''.join([
    head.tobytes(),
    comps["z_index"],
    labels_binary,
    comps["crack_codes"],
  ])

def reencode(binary:bytes, markov_model_order:int):
  return fastcrackle.reencode_markov(binary, markov_model_order)

