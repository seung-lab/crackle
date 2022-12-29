from .headers import CrackleHeader
from .codec import compress, decompress, remap, labels, nbytes, contains

class CrackleArray:
  def __init__(self, binary):
    self.binary = binary

  def __len__(self):
    return len(self.binary)

  @property
  def random_access(self):
    return True

  @property
  def size(self):
    shape = self.shape
    return shape[0] * shape[1] * shape[2]

  @property
  def nbytes(self):
    return nbytes(self.binary)

  @property
  def dtype(self):
    return header(self.binary).dtype

  @property
  def shape(self):
    head = header(self.binary)
    return (head["sx"], head["sy"], head["sz"])

  def labels(self):
    return labels(self.binary)

  def remap(self, buf, mapping, preserve_missing_labels=False):
    return CrackleArray(remap(buf, mapping, preserve_missing_labels))

  def decompress(self):
    return decompress(self.binary)

  def __contains__(self, elem):
    return contains(self.binary, elem)

  def __getitem__(self, slcs):
    slices = reify_slices(slcs, *self.shape)

    if isinstance(slcs, (slice, int)):
      slcs = (slcs,)

    while len(slcs) < 3:
       slcs += (slice(None, None, None),)

    # if self.random_access_enabled:
    #   img = decompress(self.binary, z=(slices[2].start, slices[2].stop))
    #   zslc = slice(None, None, slices[2].step)
    #   if hasattr(slcs, "__getitem__") and isinstance(slcs[2], int):
    #     zslc = 0
    #   slices = (slcs[0], slcs[1], zslc)
    #   return img[slices]
    # else:
    img = decompress(self.binary)
    return img[slcs]

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

  # First three slices are x,y,z, last is channel. 
  # Handle only x,y,z here, channel seperately
  for index, slc in enumerate(slices):
    if isinstance(slc, integer_types) or isinstance(slc, floating_types):
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