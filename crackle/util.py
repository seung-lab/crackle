import os.path
import gzip

from .encoder import compress
from .decoder import decompress

def load(filelike):
  """Load an image from a file-like object or file path."""
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
  return decompress(binary)

def save(labels, filelike):
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

