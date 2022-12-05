import numpy as np

from .lib import compute_byte_width

class FormatError(Exception):
	pass

class CrackleHeader:
  MAGIC = b'crkl'
  FORMAT_VERSION = 0
  HEADER_BYTES = 17

  def __init__(
    self, 
    data_width:int, stored_data_width:int,
    sx:int, sy:int, sz:int,
    num_label_bytes:int
  ):
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

    return CrackleHeader(
    	data_width=buffer[5],
    	stored_data_width=buffer[6],
    	sx=int.from_bytes(buffer[7:9], byteorder='little', signed=False),
    	sy=int.from_bytes(buffer[9:11], byteorder='little', signed=False),
    	sz=int.from_bytes(buffer[11:13], byteorder='little', signed=False),
    	num_label_bytes=int.from_bytes(buffer[13:17], byteorder='little', signed=False),
    )

  def tobytes(self) -> bytes:
    return b''.join([
      self.MAGIC,
      self.FORMAT_VERSION.to_bytes(1, 'little'),
      self.data_width.to_bytes(1, 'little'),
      self.stored_data_width.to_bytes(1, 'little'),
      self.sx.to_bytes(2, 'little'),
      self.sy.to_bytes(2, 'little'),
      self.sz.to_bytes(2, 'little'),
      self.num_label_bytes.to_bytes(4, 'little'),
    ])

  def index_width(self) -> int: 
    return compute_byte_width(self.sx * self.sy * self.sz)

  def depth_width(self) -> int:
    return compute_byte_width(self.sx * self.sy)

  def z_index_width(self) -> int:
    return compute_byte_width(self.sx * self.sy * 2)
