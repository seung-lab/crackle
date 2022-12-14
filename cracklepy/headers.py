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

class CrackleHeader:
  MAGIC = b'crkl'
  FORMAT_VERSION = 0
  HEADER_BYTES = 22

  def __init__(
    self, 
    label_format:int,
    crack_format:int,
    data_width:int, stored_data_width:int,
    sx:int, sy:int, sz:int,
    num_label_bytes:int
  ):
    self.label_format = label_format
    self.crack_format = crack_format
    self.data_width = int(data_width)
    self.stored_data_width = int(stored_data_width)
    self.sx = int(sx)
    self.sy = int(sy)
    self.sz = int(sz)
    self.num_label_bytes = num_label_bytes
    # should we have a field that is y/n pins?

  @classmethod
  def frombytes(kls, buffer:bytes):
    if len(buffer) < CrackleHeader.HEADER_BYTES:
    	raise FormatError(f"Bytestream too short. Got: {buffer}")
    if buffer[:4] != CrackleHeader.MAGIC:
    	raise FormatError(f"Incorrect magic number. Got: {buffer[:4]} Expected: {CrackleHeader.MAGIC}")
    if buffer[4] != CrackleHeader.FORMAT_VERSION:
    	raise FormatError(f"Wrong format version. Got: {buffer[4]} Expected: {CrackleHeader.FORMAT_VERSION}")

    values = unpack_bits(int(buffer[5]), [
      2, 2, 1, 2
    ])

    return CrackleHeader(
      label_format=values[3],
      crack_format=values[2],
    	data_width=(2 ** values[0]),
    	stored_data_width=(2 ** values[1]),
    	sx=int.from_bytes(buffer[6:10], byteorder='little', signed=False),
    	sy=int.from_bytes(buffer[10:14], byteorder='little', signed=False),
    	sz=int.from_bytes(buffer[14:18], byteorder='little', signed=False),
    	num_label_bytes=int.from_bytes(buffer[18:22], byteorder='little', signed=False),
    )

  def tobytes(self) -> bytes:
    fmt_byte = pack_bits([
      (int(np.log2(self.data_width)), 2),
      (int(np.log2(self.stored_data_width)), 2),
      (self.crack_format, 1),
      (self.label_format, 2),
    ])

    return b''.join([
      self.MAGIC,
      self.FORMAT_VERSION.to_bytes(1, 'little'),
      fmt_byte.to_bytes(1, 'little'),
      self.sx.to_bytes(4, 'little'),
      self.sy.to_bytes(4, 'little'),
      self.sz.to_bytes(4, 'little'),
      self.num_label_bytes.to_bytes(4, 'little'),
    ])

  @property
  def stored_dtype(self):
    return width2dtype[self.stored_data_width]

  @property
  def dtype(self):
    return width2dtype[self.data_width]

  def index_width(self) -> int: 
    return compute_byte_width(self.sx * self.sy * self.sz)

  def depth_width(self) -> int:
    return compute_byte_width(self.sx * self.sy)

  def z_index_width(self) -> int:
    return compute_byte_width(self.sx * self.sy * 2)

  def __repr__(self):
    return str(self.__dict__)
