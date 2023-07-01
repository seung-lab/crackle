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

def load(filelike) -> np.ndarray:
  """Load an image from a file-like object or file path."""
  return decompress(_load(filelike))

def load_numpy(filelike):
  f = io.BytesIO(_load(filelike))
  return np.load(f)

def save(labels:np.ndarray, filelike, **kwargs):
  """Save labels into the file-like object or file path."""
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

