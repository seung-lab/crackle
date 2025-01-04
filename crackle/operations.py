from typing import List, Optional, Tuple, Sequence, Union
from collections import namedtuple

import numpy as np
import fastremap
import fastcrackle

from .codec import (
  compress, decompress, labels, 
  header, raw_labels, decode_flat_labels,
  num_labels, crack_codes, components,
  reencode, background_color,
)
from .headers import CrackleHeader, CrackFormat, LabelFormat, FormatError
from .lib import width2dtype, compute_byte_width, compute_dtype

def min(binary:bytes) -> int:
  """Returns the minimum label of the crackle binary."""
  head = CrackleHeader.frombytes(binary)

  if not head.is_sorted:
    return int(np.min(labels(binary)))

  off = head.HEADER_BYTES + head.sz * 4

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

  loff = head.HEADER_BYTES + head.sz * 4

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

def remap(binary:bytes, mapping:dict, preserve_missing_labels:bool = False) -> bytes:
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
  if head.label_format == LabelFormat.PINS_VARIABLE_WIDTH:
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
    binary[CrackleHeader.HEADER_BYTES:] 
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
      binary[CrackleHeader.HEADER_BYTES:]
    ])

  return (binary, mapping)

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
    if head.label_format != LabelFormat.FLAT:
      raise ValueError("Only the FLAT label format is compatible (for now).")
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

  uniq_map = {
    u: i
    for i, u in enumerate(uniq)
  }

  component_index = []
  all_keys = []
  for binary in binaries:
    if binary is None:
      continue

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

  def synth(head, zindex, local_label_idx, keys, cracks):
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

    return b''.join([
      head.tobytes(),
      zindex.tobytes(),
      labels_binary,
      *cracks
    ])

  def synth_z_range(z_start, z_end):
    return synth(
      head, 
      all_zindex[z_start:z_end], 
      label_idx[z_start:z_end], 
      keys[label_idx_offsets[z_start]:label_idx_offsets[z_end]],
      cracks[z_start:z_end]
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
    binary[CrackleHeader.HEADER_BYTES:],
  ])

def ascontiguousarray(binary:bytes) -> bytes:
  """Convert a crackle binary to C (row-major) order."""
  head = header(binary)
  if not head.fortran_order:
    return binary

  head.fortran_order = False

  return b''.join([
    head.tobytes(),
    binary[CrackleHeader.HEADER_BYTES:],
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

  return b''.join([
    head.tobytes(),
    np.full([head.sz], len(empty_slice_crack_code), dtype=np.uint32),
    labels_binary,
    empty_slice_crack_code * head.sz,
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
  parts = decode_flat_labels(binary, head.stored_dtype, head.dtype, head.sz)
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

  return b''.join([
    head.tobytes(),
    full_parts["z_index"],
    labels_binary,
    full_parts["crack_codes"],
  ])

def add_scalar(binary:bytes, scalar:int) -> bytes:
  return operator(binary, lambda uniq: uniq + scalar)

def subtract_scalar(binary:bytes, scalar:int) -> bytes:
  return operator(binary, lambda uniq: uniq - scalar)

def multiply_scalar(binary:bytes, scalar:int) -> bytes:
  return operator(binary, lambda uniq: uniq * scalar)

def floordiv_scalar(binary:bytes, scalar:int) -> bytes:
  return operator(binary, lambda uniq: uniq // scalar)

def truediv_scalar(binary:bytes, scalar:int) -> bytes:
  return operator(binary, lambda uniq: uniq / scalar)
