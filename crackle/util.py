from typing import Optional, Union

import io
import mmap
import os
import gzip
import lzma

from .array import CrackleArray, CrackleRemoteArray
from .codec import compress, decompress
from .headers import CrackleHeader

import numpy as np

def _load(filelike, size:int = -1, allow_mmap:bool = False):
  if hasattr(filelike, 'read'):
    binary = filelike.read(size)
  elif (
    isinstance(filelike, str) 
    and os.path.splitext(filelike)[1] == '.gz'
  ):
    with gzip.open(filelike, 'rb') as f:
      binary = f.read(size)
  elif (
    isinstance(filelike, str) 
    and os.path.splitext(filelike)[1] in ('.lzma', '.xz')
  ):
    with lzma.open(filelike, 'rb') as f:
      binary = f.read(size)
  else:
    with open(filelike, 'rb') as f:
      if allow_mmap:
        binary = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
      else:
        binary = f.read(size)
  
  return binary

def load_header(filelike, **kwargs):
  """Load the header using minimal or near minimal data loading."""
  binary = _load(filelike, CrackleHeader.HEADER_BYTES)
  return CrackleHeader.frombytes(binary, **kwargs)

def load_num_labels(filelike, **kwargs):
  """Load the number of labels using near minimal data reads."""
  startpos = 0
  if hasattr(filelike, "tell"):
    startpos = filelike.tell()

  head = load_header(filelike, ignore_crc_check=kwargs.get("ignore_crc_check", False))
  readlen = head.header_bytes + head.grid_index_bytes + 16
  if hasattr(filelike, "seek"):
    filelike.seek(startpos)

  binary = _load(filelike, readlen)
  arr = CrackleArray(binary)
  return arr.num_labels()

def rload(filelike, **kwargs):
  """Load the array using a memory efficient remote interface."""
  return CrackleRemoteArray(filelike, **kwargs)

def aload(filelike, allow_mmap=False) -> CrackleArray:
  """Load a CrackleArray from a file."""
  return CrackleArray(_load(filelike, allow_mmap=allow_mmap))

def bload(filelike, allow_mmap=False) -> bytes:
  """Load the binary file."""
  return _load(filelike, allow_mmap=allow_mmap)

def load(filelike, label:Optional[int] = None) -> np.ndarray:
  """Load an image from a file-like object or file path."""
  return decompress(_load(filelike), label=label)

def load_numpy(filelike):
  f = io.BytesIO(_load(filelike))
  return np.load(f)

def save_numpy(
  arr:Union[np.ndarray, CrackleArray, bytes], 
  filelike, 
  block_size=int(200e6),
):
  if isinstance(arr, np.ndarray):
    np.save(filelike, arr)
    return

  if isinstance(arr, bytes):
    arr = CrackleArray(arr)

  if (
    isinstance(filelike, str) 
    and os.path.splitext(filelike)[1] == '.gz'
  ):
    f = gzip.open(filelike, 'wb')
  elif (
    isinstance(filelike, str) 
    and os.path.splitext(filelike)[1] in ('.lzma', '.xz')
  ):
    f = lzma.open(filelike, 'wb')
  elif isinstance(filelike, str):
    f = open(filelike, 'wb')

  head = arr.header()
  data_width = head.data_width

  np.lib.format.write_array_header_2_0(f, {
    "descr": f"<u{data_width}",
    "fortran_order": head.fortran_order,
    "shape": arr.shape,
  })

  blocks = max(int(np.ceil(arr.nbytes / block_size)), 1)
  sz = arr.shape[2]

  sz_blocks = max(int(np.ceil(sz / blocks)), 1)
  num_z_blocks = max(int(np.ceil(sz / sz_blocks)), 1)

  order = "F" if head.fortran_order else "C"

  try:
    for z_block in range(num_z_blocks):
      start = z_block * sz_blocks
      end = min((z_block+1) * sz_blocks, arr.shape[2])

      subarr = arr[:,:,start:end]
      f.write(subarr.tobytes(order))
  finally:
    if isinstance(filelike, str):
      f.close()

def save(
  labels:Union[np.ndarray, CrackleArray], 
  filelike, 
  **kwargs
):
  """Save labels into the file-like object or file path."""
  if isinstance(labels, CrackleArray):
    binary = labels.binary
  
  if (
    isinstance(filelike, str)
    and (
      filelike.endswith(".npy")
      or filelike.endswith(".npy.gz")
      or filelike.endswith(".npy.xz")
      or filelike.endswith(".npy.lzma")
    )
  ):
    return save_numpy(binary, filelike)

  if isinstance(labels, np.ndarray):
    binary = compress(labels, **kwargs)

  if hasattr(filelike, 'write'):
    filelike.write(binary)
  elif (
    isinstance(filelike, str) 
    and os.path.splitext(filelike)[1] == '.gz'
  ):
    with gzip.open(filelike, 'wb') as f:
      f.write(binary)
  elif (
    isinstance(filelike, str) 
    and os.path.splitext(filelike)[1] in ('.lzma', '.xz')
  ):
    with lzma.open(filelike, 'wb') as f:
      f.write(binary)
  else:
    with open(filelike, 'wb') as f:
      f.write(binary)

