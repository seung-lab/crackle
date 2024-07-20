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

def min(binary:bytes) -> int:
  header = CrackleHeader.frombytes(binary)

  if not head.is_sorted:
    return int(np.min(labels(binary)))

  off = header.HEADER_BYTES + header.sz * 4

  if header.label_format == LabelFormat.FLAT:
    return int.from_bytes(binary[off+8:off+8+header.stored_data_width], byteorder='little')
  else:
    bgcolor = background_color(binary)
    sdw = header.stored_data_width
    off += sdw+8
    arrmin = int.from_bytes(binary[off:off+header.stored_data_width], byteorder='little')
    if bgcolor < arrmin:
      return bgcolor
    return arrmin

def max(binary:bytes) -> int:
  header = CrackleHeader.frombytes(binary)

  if not head.is_sorted:
    return int(np.max(labels(binary)))

  loff = header.HEADER_BYTES + header.sz * 4

  if header.label_format == LabelFormat.FLAT:
    N = num_labels(binary)
    offset = loff + 8 + (N-1) * header.stored_data_width
    return int.from_bytes(binary[offset:offset+header.stored_data_width], byteorder='little')
  else:
    bgcolor = background_color(binary)
    sdw = header.stored_data_width
    N = num_labels(binary) - 1
    offset = loff + sdw + 8 + (N-1) * sdw
    arrmin = int.from_bytes(binary[offset:offset+sdw], byteorder='little')
    if bgcolor > arrmin:
      return bgcolor
    return arrmin

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
  is_sorted = np.all(all_labels[:-1] <= all_labels[1:])
  if is_sorted != head.is_sorted:
    head.is_sorted = is_sorted
    binary[:hb] = head.tobytes()

  binary[offset:offset+uniq_bytes] = list(all_labels.view(np.uint8))
  return bytes(binary)

def refit(binary:bytes) -> bytes:
  """
  Change the rendered dtype to the smallest
  dtype needed to render the image without
  loss of precision.
  """
  head = header(binary)
  dtype = fastremap.fit_dtype(head.dtype, num_labels(binary))
  head.data_width = np.dtype(dtype).itemsize
  return b''.join([ 
    head.tobytes(), 
    binary[CrackleHeader.HEADER_BYTES:] 
  ])

def renumber(binary:bytes, start=0) -> Tuple[bytes, dict]:
  """
  Renumber the array and resize the data type
  to be the smallest one to fit without loss of
  precision.
  """
  head = header(binary)
  uniq = fastremap.unique(labels(binary))
  mapping = { u: start+i for i,u in enumerate(uniq) }
  binary = refit(remap(binary, mapping))

  if not head.is_sorted:
    head.is_sorted = True
    binary[:CrackleHeader.HEADER_BYTES] = head.tobytes()

  return (binary, mapping)

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
  z_index = np.cumsum(z_index) - z_index[0]
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
  components_per_grid = np.cumsum(components_per_grid)

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

def zstack(images:Sequence[Union[np.ndarray, bytes]]) -> bytes:
  """
  Given a set of arrays or crackle compressed binaries that
  represent images of equal height and width and compatible 
  dtype, create a new crackle binary that represents a single 
  image with the input images or binaries stacked in the z-dimension.

  For example, if the input images are z-depth 3, 1, 5, 3,
  produce a new image binary that has z-size 12.

  This currently only works for flat label formats 
  and markov order 0.

  Why use this? You can iteratively build a huge array within
  limited ram.
  """
  binaries = []

  first_head = None
  sz = 0

  for binary in images:
    if isinstance(binary, np.ndarray):
      binary = compress(binary)

    head = header(binary)
    if first_head is None:
      first_head = head 

    if first_head.sx != head.sx or first_head.sy != head.sy:
      raise ValueError(
        f"All images must have the same width and height. "
        f"Expected sx={first_head.sx} sy={first_head.sy} ; Got: sx={head.sx} sy={head.sy}"
      )
    if head.label_format != LabelFormat.FLAT:
      raise ValueError("Only the FLAT label format is compatible (for now).")
    if head.markov_model_order != 0:
      raise ValueError("Markov chain encoding not currently supported.")
    if head.grid_size != first_head.grid_size:
      raise ValueError("Grid sizes must match.")
    if head.crack_format != first_head.crack_format:
      raise ValueError("All crack formats must match.")
    if head.fortran_order != first_head.fortran_order:
      raise ValueError("All binaries must be in either Fortran or C order.")
    if head.data_width != first_head.data_width:
      raise ValueError("All binaries must be the same data width.")
    if head.stored_data_width != first_head.stored_data_width:
      raise ValueError("All binaries must be the same stored data width.")
    if head.signed != first_head.signed:
      raise ValueError("All binaries must have the same sign.")

    sz += head.sz
    binaries.append(binary)

  if len(binaries) == 1:
    return binaries[0]

  first_head.sz = sz

  uniq = []
  for binary in binaries:
    uniq.extend(labels(binary))
  uniq = fastremap.unique(uniq)

  uniq_map = {
    u: i
    for i, u in enumerate(uniq)
  }

  component_index = []
  all_keys = []
  for binary in binaries:
    head = CrackleHeader.frombytes(binary)
    N = num_labels(binary)
    raw = raw_labels(binary)
    idx_bytes = head.component_width() * head.sz
    offset = 8 + N * head.stored_data_width
    component_index.append(
      np.frombuffer(raw[offset:offset + idx_bytes], dtype=f"u{head.component_width()}")
    )
    offset += idx_bytes
    key_width = compute_byte_width(N)
    keys = np.frombuffer(raw[offset:], dtype=f'u{key_width}')
    local_uniq = labels(binary)
    all_keys += [ uniq_map[key] for key in local_uniq[keys] ]

  first_head.stored_data_width = compute_byte_width(uniq.max())
  key_width = compute_byte_width(len(uniq))

  labels_binary = b''.join([
    len(uniq).to_bytes(8, 'little'),
    uniq.astype(first_head.stored_dtype).tobytes(),
    np.concatenate(component_index).tobytes(),
    np.array(all_keys, dtype=f'u{key_width}').tobytes(),
  ])

  crack_codes_lst = []
  zindex = np.zeros((first_head.sz,), dtype=np.uint32)
  z = 0
  for binary in binaries:
    for cc in crack_codes(binary):
      zindex[z] = len(cc)
      crack_codes_lst.append(cc)
      z += 1

  del binaries

  crack_binary = b''.join(crack_codes_lst)
  del crack_codes_lst

  first_head.num_label_bytes = len(labels_binary)

  return b''.join([ 
    first_head.tobytes(),
    zindex.tobytes(),
    labels_binary,
    crack_binary
  ])
