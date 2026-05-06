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

def normalize_file_ext(filename):
  filename, ext = os.path.splitext(filename)

  two_pass = ('.ckl', '.cpso')

  if ext in two_pass:
    return ext

  while True:
    filename, ext2 = os.path.splitext(filename)
    if ext2 in two_pass:
      return ext2
    elif ext2 == '':
      return ext
    ext = ext2

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

def load(filelike, label:Optional[int] = None, parallel:int = 0) -> np.ndarray:
  """Load an image from a file-like object or file path."""
  return decompress(_load(filelike), label=label, parallel=parallel)

def load_any(filename:str) -> np.ndarray:
  ext = normalize_file_ext(filename)

  if ext == ".ckl":
    image = aload(filename)
  elif ext == ".npy":
    image = load_numpy(filename)
  elif ext == ".nrrd":
    import nrrd
    image, header = nrrd.read(filename)
    if image.shape[0] == 3 and image.ndim == 3:
      image = image[...,np.newaxis]
      image = np.transpose(image, axes=[1,2,3,0])
    return image
  elif ext == ".nii":
    import nibabel as nib
    image = nib.load(filename)
    image = np.array(image.dataobj)
  elif ext in (".tif", ".tiff"):
    import tifffile
    image = tifffile.imread(filename)
  elif ext == ".cpso":
    import compresso
    image = compresso.load(filename)
  else:
    raise ValueError("Data type not supported: " + ext)

  return np.asfortranarray(image)

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

def save_nii(arr:Union[np.ndarray, bytes], path:str, affine=None):
  """Save numpy array as NIfTI. Path should end in .nii or .nii.gz"""
  import nibabel as nib
  if affine is None:
    affine = np.eye(4)

  if isinstance(arr, bytes):
    arr = crackle.decompress(arr)
  elif isinstance(arr, CrackleArray):
    arr = arr.decompress()

  img = nib.Nifti1Image(arr, affine)
  # nibabel infers gzip from .nii.gz extension automatically
  nib.save(img, path)

def save_nrrd(arr:Union[np.ndarray, bytes], path:str, compress:str = "raw"):
  """Save numpy array as NRRD."""
  import nrrd

  if isinstance(arr, bytes):
    arr = crackle.decompress(arr)
  elif isinstance(arr, CrackleArray):
    arr = arr.decompress()

  options = {}
  if compress == "gzip":
    options['encoding'] = 'gzip'
  elif compress == "bzip2":
    options['encoding'] = 'bz2'
  else:
    options['encoding'] = 'raw'

  nrrd.write(path, arr, options)

def save_tiff(arr:Union[np.ndarray, bytes], path:str, compression='zlib'):
  """
  Save numpy array as TIFF.
  compression: 'zlib' (gzip), 'lzma' (xz), or None
  """
  import tifffile

  if isinstance(arr, bytes):
    arr = crackle.decompress(arr)
  elif isinstance(arr, CrackleArray):
    arr = arr.decompress()

  tifffile.imwrite(path, arr, compression=compression)

def save_compresso(arr:Union[np.ndarray, bytes], path:str):
  import compresso

  if isinstance(arr, bytes):
    arr = crackle.decompress(arr)
  elif isinstance(arr, CrackleArray):
    arr = arr.decompress()

  compresso.save(arr, path)

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
  elif (
    isinstance(filelike, str)
    and filelike.endswith(".nrrd")
  ):
    return save_nrrd(binary, filelike)
  elif (
    isinstance(filelike, str)
    and (filelike.endswith(".tiff") or filelike.endswith(".tif"))
  ):
    return save_tiff(binary, filelike)
  elif (
    isinstance(filelike, str)
    and filelike.endswith(".cpso")
  ):
    return save_compresso(binary, filelike)

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

