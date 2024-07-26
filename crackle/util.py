from typing import Optional, Union

import io
import os
import gzip
import lzma

from .array import CrackleArray
from .codec import compress, decompress

import numpy as np

def _load(filelike):
  if hasattr(filelike, 'read'):
    binary = filelike.read()
  elif (
    isinstance(filelike, str) 
    and os.path.splitext(filelike)[1] == '.gz'
  ):
    with gzip.open(filelike, 'rb') as f:
      binary = f.read()
  elif (
    isinstance(filelike, str) 
    and os.path.splitext(filelike)[1] in ('.lzma', '.xz')
  ):
    with lzma.open(filelike, 'rb') as f:
      binary = f.read()
  else:
    with open(filelike, 'rb') as f:
      binary = f.read()
  
  return binary

def aload(filelike) -> CrackleArray:
  """Load a CrackleArray from a file."""
  return CrackleArray(_load(filelike))

def load(filelike, label:Optional[int] = None) -> np.ndarray:
  """Load an image from a file-like object or file path."""
  return decompress(_load(filelike), label=label)

def load_numpy(filelike):
  f = io.BytesIO(_load(filelike))
  return np.load(f)

def save_numpy(
  arr:Union[CrackleArray, bytes], 
  filelike, 
  block_size=int(200e6),
):
  if isinstance(filelike, str):
    f = open(filelike, "wb")

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

  order = "F" if head.fortran_order else "C"

  for z in range(blocks):
    start = z * sz_blocks
    end = min((z+1) * sz_blocks, arr.shape[2])
    
    subarr = arr[:,:,start:end]
    f.write(subarr.tobytes(order))

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
  else:
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

