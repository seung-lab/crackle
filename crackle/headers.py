from typing import List, Optional
from enum import IntEnum

import numpy as np

from .lib import (
  compute_byte_width, width2dtype,
  pack_bits, unpack_bits
)

class FormatError(Exception):
	pass

class LabelFormat(IntEnum):
  FLAT = 0
  PINS_FIXED_WIDTH = 1
  PINS_VARIABLE_WIDTH = 2

class CrackFormat(IntEnum):
  IMPERMISSIBLE = 0
  PERMISSIBLE = 1

def crc8(data:List[int]) -> int:
  # use implicit +1 representation for right shift, LSB first
  # use explicit +1 representation for left shit, MSB first
  polynomial = 0xe7 # implicit
  crc = 0xFF # detects zeroed data better than 0x0000
  for i in range(len(data)):
    crc ^= data[i]
    for k in range(8):
      if crc & 1:
        crc = (crc >> 1) ^ polynomial
      else:
        crc = crc >> 1

  return int(crc & 0xFF)


class CrackleHeader:
  MAGIC = b'crkl'
  FORMAT_VERSION = 1
  HEADER_BYTES = 29
  HEADER_BYTES_V0 = 24
  HEADER_BYTES_V1 = 29

  def __init__(
    self, 
    label_format:int,
    crack_format:int,
    data_width:int, stored_data_width:int,
    sx:int, sy:int, sz:int,
    num_label_bytes:int,
    fortran_order:bool,
    grid_size:int,
    signed:bool,
    markov_model_order:int,
    is_sorted:bool,
    format_version:int = 1,
    crc:Optional[int] = None,
  ):
    self.label_format = label_format
    self.crack_format = crack_format
    self.data_width = int(data_width)
    self.stored_data_width = int(stored_data_width)
    self.sx = int(sx)
    self.sy = int(sy)
    self.sz = int(sz)
    self.num_label_bytes = num_label_bytes
    self.fortran_order = fortran_order
    self.grid_size = int(grid_size)
    self.signed = bool(signed)
    self.markov_model_order = int(markov_model_order)
    self.is_sorted = bool(is_sorted)
    self.format_version = format_version
    self.crc = crc

  @classmethod
  def frombytes(kls, buffer:bytes, ignore_crc_check:bool = False):
    if len(buffer) < CrackleHeader.HEADER_BYTES:
      raise FormatError(f"Bytestream too short. Got: {buffer}")
    if buffer[:4] != CrackleHeader.MAGIC:
      raise FormatError(f"Incorrect magic number. Got: {buffer[:4]} Expected: {CrackleHeader.MAGIC}")

    format_version = buffer[4]
    if format_version not in [0,1]:
      raise FormatError(f"Wrong format version. Got: {format_version} Expected: <{CrackleHeader.FORMAT_VERSION}")

    values = unpack_bits(
      int.from_bytes(buffer[5:7], byteorder='little', signed=False), 
    [
      2, 2, 1, 2, 1,
      1, 4, 1
    ])

    nlabel_width = 8
    if format_version == 0:
      nlabel_width = 4
      stored_crc = None
    else:
      stored_crc = int.from_bytes(buffer[28:29], 'little')
      computed_crc = crc8(buffer[5:28])
      if not ignore_crc_check and stored_crc != computed_crc:
        raise FormatError(
          f"The header appears to be corrupted. CRC check failed. "
          f"Computed: {computed_crc} Stored: {stored_crc}"
        )

    return CrackleHeader(
      label_format=values[3],
      crack_format=values[2],
    	data_width=(2 ** values[0]),
    	stored_data_width=(2 ** values[1]),
    	sx=int.from_bytes(buffer[7:11], byteorder='little', signed=False),
    	sy=int.from_bytes(buffer[11:15], byteorder='little', signed=False),
    	sz=int.from_bytes(buffer[15:19], byteorder='little', signed=False),
      grid_size=(2 ** int(buffer[19])),
    	num_label_bytes=int.from_bytes(buffer[20:20+nlabel_width], byteorder='little', signed=False),
      fortran_order=bool(values[4]),
      signed=bool(values[5]),
      markov_model_order=int(values[6]),
      is_sorted=(not bool(values[7])),
      format_version=format_version,
      crc=stored_crc,
    )

  @property
  def header_bytes(self):
    if self.format_version == 0:
      return CrackleHeader.HEADER_BYTES_V0
    else:
      return CrackleHeader.HEADER_BYTES_V1

  @property
  def grid_index_bytes(self):
    if self.format_version == 0:
      return 4 * self.sz
    else:
      return 4 * (self.sz + 1) # plus crc32c

  def tobytes(self) -> bytes:
    fmt_byte = pack_bits([
      (int(np.log2(self.data_width)), 2),
      (int(np.log2(self.stored_data_width)), 2),
      (self.crack_format, 1),
      (self.label_format, 2),
      (self.fortran_order, 1),
    ])

    fmt_byte2 = pack_bits([
      (int(self.signed), 1),
      (int(self.markov_model_order), 4),
      (int(not self.is_sorted), 1),
    ])

    log_grid_size = int(np.log2(self.grid_size))

    fmt_ver = self.format_version
    if fmt_ver == 0 and self.num_label_bytes > np.iinfo(np.uint32).max:
      fmt_ver = 1

    if fmt_ver == 0:
      label_bytes_width = 4
    else:
      label_bytes_width = 8

    interpretable_data = b''.join([
      fmt_byte.to_bytes(1, 'little'),
      fmt_byte2.to_bytes(1, 'little'),
      self.sx.to_bytes(4, 'little'),
      self.sy.to_bytes(4, 'little'),
      self.sz.to_bytes(4, 'little'),
      log_grid_size.to_bytes(1, 'little'),
      self.num_label_bytes.to_bytes(label_bytes_width, 'little'),
    ])

    if fmt_ver == 0:
      label_bytes_width = 4
      crc = b''
    else:
      label_bytes_width = 8
      crc = crc8(interpretable_data).to_bytes(1, 'little')

    return b''.join([
      self.MAGIC,
      fmt_ver.to_bytes(1, 'little'),
      interpretable_data,
      crc
    ])

  def pin_index_width(self):
    return compute_byte_width(self.sx * self.sy * self.sz)

  def component_width(self): 
    """The size of the flat encoding type components per a grid."""
    return compute_byte_width(self.sx * self.sy)

  def num_grids(self) -> int:
    gsize = min(self.grid_size, max(self.sx, self.sy))
    ngrids = ((self.sx + gsize - 1) // gsize) * ((self.sy + gsize - 1) // gsize)
    ngrids = max(ngrids, 1)
    ngrids *= self.sz
    return int(ngrids)

  @property
  def stored_dtype(self):
    return width2dtype[self.stored_data_width]

  @property
  def dtype(self):
    return width2dtype[self.data_width]

  def index_width(self) -> int: 
    return compute_byte_width(self.sx * self.sy * self.sz)

  def depth_width(self) -> int:
    return compute_byte_width(self.sz - 1)

  def z_index_width(self) -> int:
    return 4

  @property
  def nbytes(self) -> int:
    return self.voxels() * self.data_width

  def voxels(self) -> int:
    return self.sx * self.sy * self.sz

  def compute_crc(self) -> int:
    return int.from_bytes(self.tobytes()[-2:], 'little')

  @property
  def num_markov_model_bytes(self) -> int:
    if self.markov_model_order == 0:
      return 0
    return (4 ** self.markov_model_order) * 5 // 8

  def details(self) -> str:
    label_fmt = 'FLAT'
    if self.label_format == LabelFormat.PINS_FIXED_WIDTH:
      label_fmt = 'FIXED_PINS'
    elif self.label_format == LabelFormat.PINS_VARIABLE_WIDTH:
      label_fmt = 'CONDENSED_PINS'

    return f"""
    magic:         {CrackleHeader.MAGIC}
    version:       {self.format_version}
    label fmt:     {label_fmt}
    crack fmt:     {'PERMISSIBLE' if self.crack_format == CrackFormat.PERMISSIBLE else 'IMPERMISSIBLE' }
    data width:    {self.data_width}
    stored width:  {self.stored_data_width}
    sx:            {self.sx}
    sy:            {self.sy}
    sz:            {self.sz}
    label bytes:   {self.num_label_bytes}
    fortran order: {self.fortran_order}
    grid_size:     {self.grid_size}
    crc:           {self.crc}
    ---
    BOC width:     {self.index_width()}
    z index width: {self.z_index_width()}
    """

  def __repr__(self):
    return str(self.__dict__)
