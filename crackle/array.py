from typing import Optional, Union, Any, Dict

from .headers import CrackleHeader, LabelFormat
from .codec import (
  compress, decompress, decompress_range, 
  labels, nbytes, contains, 
  header, crack_codes, num_labels,
  components, point_cloud, voxel_counts,
)
from .operations import (
  astype,
  refit, 
  renumber,
  zstack, zsplit,
  remap,
  add_scalar, subtract_scalar,
  multiply_scalar, floordiv_scalar,
)
from . import codec
from .lib import crc32c

import numpy as np

class CrackleArray:
  def __init__(self, binary:bytes, parallel:int = 0):
    self.binary = binary

    head = header(self.binary)
    self.shape = (head.sx, head.sy, head.sz)
    self.parallel = parallel

  def __len__(self):
    return len(self.binary)

  def header(self, ignore_crc_check:bool = False):
    return header(self.binary, ignore_crc_check=ignore_crc_check)

  @property
  def random_access(self):
    return True

  @property
  def size(self) -> int:
    shape = self.shape
    return shape[0] * shape[1] * shape[2]

  @property
  def ndim(self) -> int:
    return sum([ dim >= 0 for dim in self.shape ])

  @property
  def nbytes(self) -> int:
    return nbytes(self.binary)

  def copy(self):
    return CrackleArray(self.binary)

  @property
  def dtype(self):
    return header(self.binary).dtype

  def labels(self):
    return labels(self.binary)

  def num_labels(self) -> int:
    return num_labels(self.binary)

  def voxel_counts(self) -> dict:
    return voxel_counts(self.binary, parallel=self.parallel)

  def min(self) -> int:
    return codec.min(self.binary)

  def max(self) -> int:
    return codec.max(self.binary)

  def remap(self, mapping:Dict[int,int], preserve_missing_labels:bool = False):
    return CrackleArray(remap(self.binary, mapping, preserve_missing_labels))

  def refit(self):
    return CrackleArray(refit(self.binary))

  def astype(self, dtype):
    return CrackleArray(astype(self.binary, dtype))

  def renumber(self, start:int = 0) -> bytes:
    return CrackleArray(renumber(self.binary, start))

  def numpy(self, *args, **kwargs) -> np.ndarray:
    return self.decompress(*args, **kwargs)

  def decompress(self, label:Optional[int] = None) -> bytes:
    return decompress(self.binary, label, parallel=self.parallel)

  def point_cloud(
    self, 
    label:Optional[int] = None,
    z_start:int = -1,
    z_end:int = -1,
  ) -> np.ndarray:
    return point_cloud(
      self.binary, label, 
      z_start=z_start, z_end=z_end, 
      parallel=self.parallel
    )

  def save(self, filelike):
    import crackle.util
    return crackle.util.save(self, filelike)

  def __eq__(self, other:Union[int, Any]) -> bool:
    if isinstance(other, int):
      return self.min() == other and self.max() == other
    elif isinstance(other, CrackleArray):
      return self.binary == other.binary
    else:
      raise TypeError(f"Type {type(other)} is not supported.")

  def __add__(self, other:int) -> bytes:
    return CrackleArray(add_scalar(self.binary, other))

  def __radd__(self, other:int) -> bytes:
    return self.__add__(other)

  def __sub__(self, other:int) -> bytes:
    return CrackleArray(subtract_scalar(self.binary, other))

  def __rsub__(self, other:int) -> bytes:
    return self.__sub__(other)

  def __mul__(self, other:int) -> bytes:
    return CrackleArray(multiply_scalar(self.binary, other))

  def __rmul__(self, other:int) -> bytes:
    return self.__mul__(other)

  def __floordiv__(self, other:int) -> bytes:
    return CrackleArray(floordiv_scalar(self.binary, other))

  def __rfloordiv__(self, other:int) -> bytes:
    return self.__floordiv__(other)

  def __contains__(self, elem:int) -> bool:
    return contains(self.binary, elem)

  def __getitem__(self, slcs) -> np.ndarray:
    if slcs == (Ellipsis, np.newaxis):
      self.shape = self.shape + (1,)
      return self

    slices = reify_slices(slcs, *self.shape[:3])

    if isinstance(slcs, (slice, int)):
      slcs = (slcs,)

    while len(slcs) < 3:
       slcs += (slice(None, None, None),)

    img = decompress_range(
      self.binary, 
      slices[2].start, slices[2].stop, 
      parallel=self.parallel
    )
    zslc = slice(None, None, slices[2].step)
    if isinstance(slcs[2], (int, np.integer)):
      zslc = 0

    slices = (slcs[0], slcs[1], zslc)
    cutout = img[slices]

    for _ in range(self.ndim - 3):
      cutout = cutout[..., np.newaxis]

    return cutout

  def __setitem__(self, slcs, data):
    if slcs == (Ellipsis, np.newaxis):
      self.shape = self.shape + (1,)
      return self

    slices = reify_slices(slcs, *self.shape[:3])

    if isinstance(slcs, (slice, int)):
      slcs = (slcs,)

    slice_all = slice(None, None, None)

    while len(slcs) < 3:
       slcs += (slice_all,)

    head = self.header()
    sz = slices[2].stop - slices[2].start

    # This could be alternatively written as a remap operation
    if isinstance(data, (int, float)):
      data = np.full(
        [ self.shape[0], self.shape[1], sz ],
        data,
        dtype=head.dtype,
        order=('F' if head.fortran_order else 'C'),
      )

    if slices[0] != slice(0, self.shape[0], 1) or slices[1] != slice(0, self.shape[1], 1):
      tmp_data = self[:,:,slices[2].start:slices[2].stop]
      tmp_data[(slices[0], slices[1])] = data
      data = tmp_data

    if data.shape[2] != sz:
      raise ValueError(f"{data.shape[2]} did not match slice dimensions.")

    data_binary = compress(data.astype(head.dtype, copy=False))
    
    if slices[2] == slice(0, self.shape[2], 1):
      self.binary = data_binary
      return

    (before_0, mid_0, after_0) = zsplit(self.binary, slices[2].start)
    (before_1, mid_1, after_1) = zsplit(self.binary, slices[2].stop)
    
    del mid_0
    del after_0
    del before_1

    self.binary = zstack([
      before_0,
      data_binary,
      mid_1,
      after_1,
    ])

class CrackleRemoteArray(CrackleArray):
  """EXPERIMENTAL DO NOT RELY ON THIS INTERFACE."""
  def __init__(self, cloudpath:str, ignore_header_crc_check:bool = False):
    from cloudfiles import CloudFile
    self.cloudpath = cloudpath
    self.cf = CloudFile(cloudpath)
    self.header_binary = self.cf[:CrackleHeader.HEADER_BYTES]
    self.header = header(self.header_binary, ignore_crc_check=ignore_header_crc_check)

    self.z_index = None
    self.labels_binary = None
    self.markov_model = None

  def labels(self):
    binary = self._synthetic_crackle_file(0, b'')
    return CrackleArray(binary).labels()

  def num_labels(self):
    hb = self.header.header_bytes
    offset = hb + self.header.sz * 4
    sdw = self.header.stored_data_width

    if self.header.label_format == LabelFormat.FLAT:
      num_labels_binary = self.cf[offset:offset+8]
    else:
      num_labels_binary = self.cf[offset+sdw:offset+sdw+8]

    return int.from_bytes(num_labels_binary, 'little')

  def __contains__(self, elem:int):
    binary = self._synthetic_crackle_file(0, b'')
    return elem in CrackleArray(binary)
  
  def fetch_z_index_labels_markov_model(self):
    hb = self.header.header_bytes
    z_offset = self.header.grid_index_bytes
    offset = (
      z_offset 
      + self.header.num_label_bytes 
      + self.header.num_markov_model_bytes
    )
    binary = self.cf[hb:hb+offset]

    z_index = np.frombuffer(binary[:z_offset], dtype=np.uint32)
    lo = z_offset+self.header.num_label_bytes
    labels_binary = binary[z_offset:lo]
    markov_binary = binary[lo:lo+self.header.num_markov_model_bytes]
    del binary

    z_index = np.cumsum(z_index)
    z_index = np.concatenate(([ 0 ], z_index))
    z_index += (
      hb
      + self.header.num_label_bytes 
      + self.header.grid_index_bytes
      + self.header.num_markov_model_bytes
    )
    z_index = z_index.astype(np.uint64, copy=False)
    return (z_index, labels_binary, markov_binary)

  def fetch_markov_model(self):
    if self.header.markov_model_order == 0:
      return b''

    hb = self.header.header_bytes
    off = hb + self.header.grid_index_bytes
    off += self.header.num_label_bytes
    return self.cf[off:off+self.header.num_markov_model_bytes]

  def fetch_all_labels(self) -> bytes:
    hb = self.header.header_bytes
    off = hb + self.header.grid_index_bytes
    return self.cf[off:off+self.header.num_label_bytes]

  def fetch_crack_code(self, z:int) -> bytes:
    return self.cf[self.z_index[z]:self.z_index[z+1]]

  def _synthetic_crackle_file(self, z:int, crackcode:bytes, labels_binary:Optional[bytes] = None) -> bytes:
    zindex = np.zeros((self.header.sz,), dtype=np.uint32)
    zindex[z] = len(crackcode)

    if labels_binary is None:
      labels_binary = self.labels_binary

    grid_index = zindex.tobytes()
    if self.header.format_version > 0:
      grid_index += crc32c(grid_index).to_bytes(4, 'little')

    return b''.join([ 
      self.header_binary,
      grid_index,
      labels_binary,
      self.markov_model,
      crackcode
    ])

  def __getitem__(self, z:int) -> np.ndarray:
    if self.z_index is None:
      (
        self.z_index, 
        self.labels_binary,
        self.markov_model
      ) = self.fetch_z_index_labels_markov_model()

    crackcode = self.fetch_crack_code(z)
    binary = self._synthetic_crackle_file(z, crackcode)
    return CrackleArray(binary)[:,:,z]

def reify_slices(slices, sx, sy, sz):
  """
  Convert free attributes of a slice object 
  (e.g. None (arr[:]) or Ellipsis (arr[..., 0]))
  into bound variables in the context of this
  bounding box.

  That is, for a ':' slice, slice.start will be set
  to the value of the respective minpt index of 
  this bounding box while slice.stop will be set 
  to the value of the respective maxpt index.

  Example:
    reify_slices( (np._s[:],) )
    
    >>> [ slice(-1,1,1), slice(-2,2,1), slice(-3,3,1) ]

  Returns: [ slice, ... ]
  """
  ndim = 3
  minpt = (0,0,0)
  maxpt = (sx,sy,sz)

  integer_types = (int, np.integer)
  floating_types = (float, np.floating)

  if isinstance(slices, integer_types) or isinstance(slices, floating_types):
    slices = [ slice(int(slices), int(slices)+1, 1) ]
  elif isinstance(slices, slice):
    slices = [ slices ]
  elif slices is Ellipsis:
    slices = []

  slices = list(slices)

  for index, slc in enumerate(slices):
    if slc is Ellipsis:
      fill = ndim - len(slices) + 1
      slices = slices[:index] +  (fill * [ slice(None, None, None) ]) + slices[index+1:]
      break

  while len(slices) < ndim:
    slices.append( slice(None, None, None) )

  while len(slices) > ndim and slices[-1] == slice(None, None, None):
    slices.pop()

  # First three slices are x,y,z, last is channel. 
  # Handle only x,y,z here, channel seperately
  for index, slc in enumerate(slices):
    if isinstance(slc, integer_types) or isinstance(slc, floating_types):
      slc = int(slc)
      if slc < 0:
        slc += maxpt[index]
      slices[index] = slice(int(slc), int(slc)+1, 1)
    elif slc == Ellipsis:
      raise ValueError("More than one Ellipsis operator used at once.")
    else:
      start = 0 if slc.start is None else slc.start
      end = maxpt[index] if slc.stop is None else slc.stop 
      step = 1 if slc.step is None else slc.step

      if step < 0:
        raise ValueError(f'Negative step sizes are not supported. Got: {step}')

      if start < 0: # this is support for negative indicies
        start = maxpt[index] + start         
      check_bounds(start, minpt[index], maxpt[index])
      if end < 0: # this is support for negative indicies
        end = maxpt[index] + end
      check_bounds(end, minpt[index], maxpt[index])

      slices[index] = slice(start, end, step)

  return slices

def clamp(val, low, high):
  return min(max(val, low), high)

def check_bounds(val, low, high):
  if val > high or val < low:
    raise ValueError(f'Value {val} cannot be outside of inclusive range {low} to {high}')
  return val