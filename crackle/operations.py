from typing import List, Optional, Tuple, Sequence, Union, Dict
from collections import namedtuple, defaultdict
import multiprocessing as mp
import warnings

import numpy as np
import fastremap
import fastcrackle

from .codec import (
  compress, decompress, decompress_range, labels, 
  header, raw_labels, decode_flat_labels,
  decode_condensed_pins, decode_condensed_pins_components,
  num_labels, crack_codes, components,
  reencode, background_color, crack_crcs, labels_crc,
)
from .headers import CrackleHeader, CrackFormat, LabelFormat, FormatError
from .lib import width2dtype, compute_byte_width, compute_dtype, crc32c

def min(binary:bytes) -> int:
  """Returns the minimum label of the crackle binary."""
  head = CrackleHeader.frombytes(binary)

  if not head.is_sorted:
    return int(np.min(labels(binary)))

  off = head.header_bytes + head.grid_index_bytes

  if head.label_format == LabelFormat.FLAT:
    return int.from_bytes(binary[off+8:off+8+head.stored_data_width], byteorder='little')
  else:
    bgcolor = background_color(binary)
    sdw = head.stored_data_width
    off += sdw+8
    arrmin = int.from_bytes(binary[off:off+head.stored_data_width], byteorder='little')
    if bgcolor < arrmin:
      return bgcolor
    return arrmin

def max(binary:bytes) -> int:
  """Returns the maximum label of the crackle binary."""
  head = CrackleHeader.frombytes(binary)

  if not head.is_sorted:
    return int(np.max(labels(binary)))

  loff = head.header_bytes + head.grid_index_bytes

  if head.label_format == LabelFormat.FLAT:
    N = num_labels(binary)
    offset = loff + 8 + (N-1) * head.stored_data_width
    return int.from_bytes(binary[offset:offset+head.stored_data_width], byteorder='little')
  else:
    bgcolor = background_color(binary)
    sdw = head.stored_data_width
    N = num_labels(binary) - 1
    offset = loff + sdw + 8 + (N-1) * sdw
    arrmin = int.from_bytes(binary[offset:offset+sdw], byteorder='little')
    if bgcolor > arrmin:
      return bgcolor
    return arrmin

def remap(
  binary:bytes, 
  mapping:dict, 
  preserve_missing_labels:bool = False,
  in_place:bool = False,
  parallel:int = 0,
) -> bytes:
  """
  Remap the labels in a crackle bystream without decompressing.
  
  binary: A Crackle bytestream
  mapping: A dict mapping old labels to new labels.
  preserve_missing_labels: By default we presume all values
    will be remapped. If a value in the binary does not have
    a corresponding mapping key, it will raise a KeyError. If
    True, just leave it be and encode all the labels that are 
    in the mapping.
  in_place: modify the bytestream in place. note: this will
    even modify "immutable" bytes objects.
  """
  if not in_place:
    binary = bytearray(binary)
  fastcrackle.remap(binary, mapping, preserve_missing_labels, parallel)
  return binary

def astype(binary:bytes, dtype) -> bytes:
  """
  Change the rendered dtype to the smallest
  dtype needed to render the image without
  loss of precision.
  """
  head = header(binary)
  head.data_width = np.dtype(dtype).itemsize
  return b''.join([ 
    head.tobytes(), 
    binary[head.header_bytes:] 
  ])

def refit(binary:bytes) -> bytes:
  """
  Change the rendered dtype to the smallest
  dtype needed to render the image without
  loss of precision.
  """
  head = header(binary)
  dtype = fastremap.fit_dtype(head.dtype, max(binary))
  return astype(binary, dtype)

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
    binary = b''.join([
      head.tobytes(),
      binary[head.header_bytes:]
    ])

  return (binary, mapping)

def _zstack_flat_labels(
  uniq:np.ndarray, binaries:List[bytes]
) -> bytes:
  """
  Convert a list of crackle binaries with flat label 
  type into a merged labels binary section.
  """
  component_index = []
  all_keys = []

  uniq_map = {
    u: i
    for i, u in enumerate(uniq)
  }

  first_head = CrackleHeader.frombytes(binaries[0])
  first_head.stored_data_width = compute_byte_width(uniq.max())
  key_width = compute_byte_width(len(uniq))

  for binary in binaries:
    if binary is None:
      continue
    head = CrackleHeader.frombytes(binary)
    elements = decode_flat_labels(head, binary)
    component_index.append(elements["components_per_grid"])
    
    local_uniq = elements["unique"]
    cc_map = elements["cc_map"]
    remap = np.array([ uniq_map[key] for key in local_uniq  ], dtype=f"u{key_width}")
    all_keys.append(remap[cc_map])
  
  # labels binary
  return b''.join([
    len(uniq).to_bytes(8, 'little'),
    uniq.astype(first_head.stored_dtype, copy=False).tobytes(),
    np.concatenate(component_index).tobytes(),
    np.concatenate(all_keys).tobytes(),
  ])

def _zstack_pins(
  uniq:np.ndarray, 
  binaries:List[bytes],
) -> bytes:
  if len(binaries) <= 1:
    return binaries

  binaries = [
    binary
    for binary in binaries 
    if binary is not None 
  ]

  first_head = CrackleHeader.frombytes(binaries[0])

  component_index = []

  # fmt: 
  # bg color, N unique, unique, cc_per_grid, fmt_byte, pins
  first_bgcolor = background_color(binaries[0])
  component_offset = 0
  z = 0
  sxy = first_head.sx * first_head.sy

  all_pins = defaultdict(list)
  all_single_labels = defaultdict(list)

  for binary in binaries:
    bgcolor = background_color(binary)
    if bgcolor != first_bgcolor:
      raise ValueError(
        f"Unable to stack pins with different background colors. "
        f"Got: {first_bgcolor} and {bgcolor}"
      )
    elems = decode_condensed_pins_components(binary)
    cpg = elems["components_per_grid"]
    component_index.append(cpg)
    del elems
    pins, single_labels = decode_condensed_pins(binary)
    for label, cc_labels in single_labels.items():
      cc_labels += int(component_offset)
      all_single_labels[label].extend(list(cc_labels))

    component_offset += int(np.sum(cpg))

    PinTuple = namedtuple('Pin', ['index', 'depth'])

    for label in pins.keys():
      all_pins[label] += [
        PinTuple(pin.index + z * sxy, pin.depth)
        for pin in pins[label]
      ]

    head = CrackleHeader.frombytes(binary)
    z += head.sz

  maxfn = __builtins__["max"] # name collision
  num_pins = maxfn([ len(v) for v in all_pins.values() ])

  max_depth = maxfn((
    pin.depth
    for label, pins in all_pins.items() 
    for pin in pins
  ))
  max_ccl = maxfn((
    ccl
    for label, ccls in all_single_labels.items()
    for ccl in ccls 
  ))

  num_pins_width = compute_byte_width(num_pins)
  depth_width = compute_byte_width(max_depth)
  cc_label_width = compute_byte_width(max_ccl)

  fmt_byte = (
    int(np.log2(num_pins_width))
    | (int(np.log2(depth_width)) << 2)
    | (int(np.log2(cc_label_width)) << 4)
  )

  # pins: | num_pins | INDEX_0 | INDEX_1 | ... | INDEX_N 
  #       | DEPTH_0 | DEPTH_1 | ... | DEPTH_N | 
  #         num_single_labels | CC 0 | CC 1 | ... | CC N |

  index_width = first_head.pin_index_width()

  pin_binaries = []
  for label in uniq:
    if label == first_bgcolor:
      continue
    pinset = all_pins[label]
    singles = all_single_labels[label]
    pinset.sort(key=lambda a: a.index)

    indices = np.array([ p.index for p in pinset ], dtype=f"u{index_width}")
    indices = np.diff(indices, prepend=0).astype(f"u{index_width}")

    depths = np.array([ p.depth for p in pinset ], dtype=f"u{depth_width}")

    single_labels = np.array(all_single_labels[label], dtype=f"u{cc_label_width}")
    single_labels.sort()
    single_labels = np.diff(single_labels, prepend=0)
    single_labels = single_labels.astype(f"u{cc_label_width}")

    pin_section = b''.join([
      len(pinset).to_bytes(num_pins_width, 'little'),
      indices.tobytes(),
      depths.tobytes(),
      len(single_labels).to_bytes(num_pins_width, 'little'),
      single_labels.tobytes(),
    ])
    pin_binaries.append(pin_section)

  uniq = uniq[uniq != first_bgcolor]

  # labels binary
  return b''.join([
    int(first_bgcolor).to_bytes(first_head.stored_data_width, 'little'),
    len(uniq).to_bytes(8, 'little'),
    uniq.astype(first_head.stored_dtype, copy=False).tobytes(),
    np.concatenate(component_index).tobytes(),
    fmt_byte.to_bytes(1, 'little'),
    *pin_binaries
  ])


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
    if binary is None:
      continue

    if isinstance(binary, np.ndarray):
      binary = compress(binary)
    else:
      binary = reencode(binary, markov_model_order=0)

    head = header(binary)
    if first_head is None:
      first_head = head 

    if first_head.fortran_order:
      binary = asfortranarray(binary)
    else:
      binary = ascontiguousarray(binary)

    if first_head.sx != head.sx or first_head.sy != head.sy:
      raise ValueError(
        f"All images must have the same width and height. "
        f"Expected sx={first_head.sx} sy={first_head.sy} ; Got: sx={head.sx} sy={head.sy}"
      )
    if first_head.label_format != head.label_format:
      raise ValueError(f"Label formats must match. First: {first_head.label_format} Got: {head.label_format}")

    if head.grid_size != first_head.grid_size:
      raise ValueError("Grid sizes must match.")
    if head.crack_format != first_head.crack_format:
      raise ValueError("All crack formats must match.")
    if head.data_width != first_head.data_width:
      raise ValueError("All binaries must be the same data width.")
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

  if first_head.label_format == LabelFormat.FLAT:
    labels_binary = _zstack_flat_labels(uniq, binaries)
  elif first_head.label_format == LabelFormat.PINS_VARIABLE_WIDTH:
    labels_binary = _zstack_pins(uniq, binaries)
  else:
    raise ValueError(f"Unsupported label format: {first_head.label_format}")

  crack_codes_lst = []
  zindex = np.zeros((sz,), dtype=np.uint32)
  z = 0
  for binary in binaries:
    for cc in crack_codes(binary):
      zindex[z] = len(cc)
      crack_codes_lst.append(cc)
      z += 1

  grid_index_binary = zindex.tobytes()
  if first_head.format_version > 0:
    computed_crc32c = crc32c(grid_index_binary)
    grid_index_binary += computed_crc32c.to_bytes(4, 'little')

  crcs_binary = b''
  if first_head.format_version > 0:
    crcs = []
    for binary in binaries:
      crcs.append(crack_crcs(binary))
    crcs_binary = np.concatenate(crcs).tobytes()
    
  del zindex
  del binaries

  crack_binary = b''.join(crack_codes_lst)
  del crack_codes_lst

  first_head.num_label_bytes = len(labels_binary)

  labels_crc = b''
  if first_head.format_version > 0:
    labels_crc = crc32c(labels_binary).to_bytes(4, 'little')

  return b''.join([ 
    first_head.tobytes(),
    grid_index_binary,
    labels_binary,
    crack_binary,
    labels_crc,
    crcs_binary
  ])

def _zsplit_helper(binary:bytes):
  head = header(binary)

  if head.label_format != LabelFormat.FLAT:
    raise ValueError("Label format not currently supported.")

  uniq = labels(binary)
  raw = raw_labels(binary)

  N = num_labels(binary)
  idx_bytes = head.component_width() * head.sz
  offset = 8 + N * head.stored_data_width
  label_idx = np.frombuffer(raw[offset:offset + idx_bytes], dtype=f"u{head.component_width()}")
  offset += idx_bytes
  key_width = compute_byte_width(N)
  keys = np.frombuffer(raw[offset:], dtype=f'u{key_width}')
  
  label_idx_offsets = np.concatenate([ [0], label_idx ])
  label_idx_offsets = np.cumsum(label_idx_offsets)

  all_zindex = np.frombuffer(components(binary)["z_index"], dtype=np.uint32)

  cracks = crack_codes(binary)
  all_crack_crcs = crack_crcs(binary)

  def synth(head, zindex, local_label_idx, keys, cracks, sub_crack_crcs):
    local_uniq = fastremap.unique(uniq[keys])
    local_uniq_map = { u: i for i, u in enumerate(local_uniq) }
    remapped_keys = [ local_uniq_map[k] for k in uniq[keys] ]

    key_width = compute_byte_width(len(local_uniq))
    head.stored_data_width = compute_byte_width(local_uniq.max())

    labels_binary = b''.join([
      len(local_uniq).to_bytes(8, 'little'),
      local_uniq.astype(head.stored_dtype).tobytes(),
      local_label_idx.tobytes(),
      np.array(remapped_keys, dtype=f'u{key_width}').tobytes(),
    ])

    head.sz = len(cracks)
    head.num_label_bytes = len(labels_binary)

    grid_index = zindex.tobytes()
    labels_crc = b''
    crack_crcs_binary = b''
    if head.format_version > 0:
      grid_index += crc32c(grid_index).to_bytes(4, 'little')
      labels_crc = crc32c(labels_binary).to_bytes(4, 'little')
      crack_crcs_binary = sub_crack_crcs.tobytes()

    return b''.join([
      head.tobytes(),
      grid_index,
      labels_binary,
      *cracks,
      labels_crc,
      crack_crcs_binary,
    ])

  def synth_z_range(z_start:int, z_end:int) -> bytes:
    tmp_crack_crcs = []
    if head.format_version > 0:
      tmp_crack_crcs = all_crack_crcs[z_start:z_end]

    return synth(
      head, 
      all_zindex[z_start:z_end], 
      label_idx[z_start:z_end], 
      keys[label_idx_offsets[z_start]:label_idx_offsets[z_end]],
      cracks[z_start:z_end],
      tmp_crack_crcs,
    )

  return synth_z_range

def zsplit(binary:bytes, z:int) -> Tuple[bytes, bytes, bytes]:
  """
  Given a crackle binary, split that binary at a given z
  into before, middle, and after binaries.
  
  Combined with zstack, this gives you a way to start editing
  binaries without full decompression.
  """
  head = header(binary)

  if z < 0 or z >= head.sz:
    raise ValueError(f"{z} is outside the range 0 to {head.sz}.")

  if head.sz == 1 and z_start == 0:
    return (b'', binary, b'')

  crt = _zsplit_helper(binary)
  before = crt(0, z)
  middle = crt(z, z+1)
  after = crt(z+1, head.sz)

  return (before, middle, after)

def zshatter(binary:bytes) -> List[bytes]:
  """
  Given a crackle binary, split that binary into
  individual z slices.
  
  Combined with zstack, this gives you a way to edit
  binaries without full decompression.
  """
  head = header(binary)
  crt = _zsplit_helper(binary)
  return [
    crt(z, z+1)
    for z in range(head.sz)
  ]

def asfortranarray(binary:bytes) -> bytes:
  """Convert a crackle binary to Fortran (column-major) order."""
  head = header(binary)
  if head.fortran_order:
    return binary

  head.fortran_order = True

  return b''.join([
    head.tobytes(),
    binary[head.header_bytes:],
  ])

def ascontiguousarray(binary:bytes) -> bytes:
  """Convert a crackle binary to C (row-major) order."""
  head = header(binary)
  if not head.fortran_order:
    return binary

  head.fortran_order = False

  return b''.join([
    head.tobytes(),
    binary[head.header_bytes:],
  ])

def full(shape, fill_value, dtype=None, order='C') -> bytes:
  """
  Create a crackle binary that represents an array
  filled with a single value. Arguments are identical
  to np.full.
  """
  if dtype is None:
    dtype = np.array(fill_value).dtype

  head = CrackleHeader(
    label_format=LabelFormat.FLAT,
    crack_format=CrackFormat.IMPERMISSIBLE,
    data_width=np.dtype(dtype).itemsize, 
    stored_data_width=compute_byte_width(fill_value),
    sx=shape[0], 
    sy=shape[1], 
    sz=shape[2],
    num_label_bytes=0,
    fortran_order=(order == 'F'),
    grid_size=int(2 ** 31),
    signed=(fill_value < 0),
    markov_model_order=0,
    is_sorted=True,
  )

  labels_binary = b''.join([
    int(1).to_bytes(8, 'little'),
    np.array([ fill_value ], dtype=head.stored_dtype).tobytes(),
    np.ones([ head.sz ], dtype=f'u{head.component_width()}'),
    np.zeros([ head.sz ], dtype=np.uint8),
  ])

  head.num_label_bytes = len(labels_binary)
  head.is_sorted = True

  empty_slice_crack_code = b'\x01\x00\x00\x00\x00'

  grid_index = np.full([head.sz], len(empty_slice_crack_code), dtype=np.uint32).tobytes()
  grid_index += crc32c(grid_index).to_bytes(4, 'little')

  labels_crc_binary = crc32c(labels_binary).to_bytes(4, 'little')
  crack_crc_single = crc32c(np.zeros(shape[:2], dtype=np.uint32))
  crack_crcs_binary = np.array([crack_crc_single] * shape[2], dtype=np.uint32).tobytes()

  return b''.join([
    head.tobytes(),
    grid_index,
    labels_binary,
    empty_slice_crack_code * head.sz,
    labels_crc_binary,
    crack_crcs_binary,
  ])

def zeros(shape, dtype=None, order="C"):
  """
  Return a new array of given shape and type, filled with zeros.
  """
  return full(shape, 0, dtype, order)

def ones(shape, dtype=None, order="C"):
  """
  Return a new array of given shape and type, filled with ones.
  """
  return full(shape, 1, dtype, order)

def operator(binary:bytes, fn) -> bytes:
  head = header(binary)
  parts = decode_flat_labels(head, binary)
  parts["unique"] = fn(parts["unique"])

  head.stored_data_width = compute_byte_width(parts["unique"][-1])

  labels_binary = b''.join([
    len(parts["unique"]).to_bytes(8, 'little'),
    parts["unique"].astype(head.stored_dtype, copy=False).tobytes(),
    parts["components_per_grid"].tobytes(),
    parts["cc_map"].tobytes(),
  ])

  full_parts = components(binary)
  head.num_label_bytes = len(labels_binary)

  labels_crc_binary = b''
  crack_crcs_binary = b''
  if head.format_version > 0:
    labels_crc_binary = crc32c(labels_binary).to_bytes(4, 'little')
    crack_crcs_binary = crack_crcs(binary).tobytes()

  return b''.join([
    head.tobytes(),
    full_parts["z_index"],
    labels_binary,
    full_parts["crack_codes"],
    labels_crc_binary,
    crack_crcs_binary,
  ])

def add_scalar(binary:bytes, scalar:int) -> bytes:
  if scalar == 0:
    return binary
  return operator(binary, lambda uniq: uniq + scalar)

def subtract_scalar(binary:bytes, scalar:int) -> bytes:
  if scalar == 0:
    return binary
  return operator(binary, lambda uniq: uniq - scalar)

def multiply_scalar(binary:bytes, scalar:int) -> bytes:
  if scalar == 1:
    return binary
  return operator(binary, lambda uniq: uniq * scalar)

def floordiv_scalar(binary:bytes, scalar:int) -> bytes:
  if scalar == 1:
    return binary
  return operator(binary, lambda uniq: uniq // scalar)

def truediv_scalar(binary:bytes, scalar:int) -> bytes:
  if scalar == 1:
    return binary
  return operator(binary, lambda uniq: uniq / scalar)

def recompress(
  binary:bytes, 
  memory_target:int = int(4e9), 
  allow_pins:bool = False,
) -> bytes:
  """
  After e.g. remapping a volume, you might want to recompress it
  to eliminate false boundaries. This function performs the recompression
  cycle relatively efficiently in memory.

  You can also use this function to toggle between pin and flat label
  encoding. However, you must have sufficient memory for at least 3-4 
  uncompressed sections (ideally more) at once for this to be "profitable".
  """
  head = header(binary)
  # approximate RAM needed to encode/decode a full z-section
  # section + CCL + VCG
  section_bytes = head.sx * head.sy * (head.data_width + 4 + 1)

  min_memory_request = 2 * len(binary) + section_bytes

  if min_memory_request >= memory_target:
    warnings.warn(
      f"Memory required is potentially more than the memory target specified. "
      f"Sometimes recompression will result in substantial additional savings "
      f"so this isn't a hard error. Recompression will require storing approximately"
      f"thirce the same binary size in RAM plus additional memory per a parallel section.\n"
      f"memory target: {memory_target} bytes\n"
      f"requested: {min_memory_request} bytes"
    )

  parallel = __builtins__["max"](memory_target - len(binary), 0) // section_bytes
  parallel = __builtins__["max"](parallel, 1)
  parallel = __builtins__["min"](parallel, mp.cpu_count())

  bgcolor = min(binary) # TODO: could be properly calculated from cc keys

  binaries = []
  for z in range(0, head.sz, parallel):
    z_end = __builtins__["min"](z+parallel, head.sz)
    labels = decompress_range(binary, z_start=z, z_end=z_end, parallel=parallel)
    binaries.append(
      compress(labels, allow_pins=allow_pins, bgcolor=bgcolor)
    )

  return zstack(binaries)


