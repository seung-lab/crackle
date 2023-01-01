import os.path
import gzip

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

def save(labels:np.ndarray, filelike):
  """Save labels into the file-like object or file path."""
  binary = compress(labels)
  if hasattr(filelike, 'write'):
    filelike.write(binary)
  elif (
    isinstance(filelike, str) 
    and os.path.splitext(filelike)[1] == '.gz'
  ):
    with gzip.open(filelike, 'wb') as f:
      f.write(binary)
  else:
    with open(filelike, 'wb') as f:
      f.write(binary)

