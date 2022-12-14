import numpy as np

class DisjointSet:
  def __init__(self):
    self.data = {} 
  def __len__(self):
    return len(self.data)
  def add(self, x):
    self.data[x] = x
    return x
  def find(self, x):
    if not x in self.data:
      return None
    i = self.data[x]
    while i != self.data[i]:
      self.data[i] = self.data[self.data[i]]
      i = self.data[i]
    return i
  def unify(self, x, y):
    i = self.find(x)
    j = self.find(y)
    if i is None:
      i = self.add(x)
    if j is None:
      j = self.add(y)

    if i < j:
      self.data[j] = i
    else:
      self.data[i] = j

def color_connectivity_graph(vcg):
  sx, sy = vcg.shape

  equivalences = DisjointSet()
  out = np.zeros((sx,sy), dtype=np.uint32)

  new_label = 0
  equivalences.add(new_label)
  for x in range(sx):
    if x > 0 and vcg[x,0] & 0b0010 == 0:
      new_label += 1
      equivalences.add(new_label)
    out[x,0] = new_label

  for y in range(1,sy):
    for x in range(sx):
      if x > 0 and (vcg[x,y] & 0b0010):
        out[x,y] = out[x-1,y]
        if y > 0 and (vcg[x,y-1] & 0b0010) == 0 and (vcg[x,y] & 0b1000) > 0:
          equivalences.unify(out[x,y], out[x,y-1])
      elif y > 0 and (vcg[x,y] & 0b1000):
        out[x,y] = out[x,y-1]
      else:
        new_label += 1
        out[x,y] = new_label
        equivalences.add(new_label)

  return relabel(out, equivalences)

def connected_components(labels):
  """4 connected CCL"""
  sx, sy = labels.shape

  equivalences = DisjointSet()
  out = np.zeros((sx,sy), dtype=np.uint32)

  new_label = -1
  for y in range(sy):
    for x in range(sx):
      if x > 0 and labels[x-1,y] == labels[x,y]:
        out[x,y] = out[x-1,y]
        if (
          y > 0 
          and (labels[x,y] != labels[x-1,y-1]) 
          and (labels[x,y] == labels[x,y-1])
        ):
          equivalences.unify(out[x,y], out[x,y-1])
      elif y > 0 and labels[x,y] == labels[x,y-1]:
        out[x,y] = out[x,y-1]
      else:
        new_label += 1
        out[x,y] = new_label
        equivalences.add(new_label)
        
  return relabel(out, equivalences)

def relabel(out, equivalences):
  sx, sy = out.shape
  nlabels = len(equivalences)

  next_label = 0
  renumber = [None] * nlabels
  for i in range(nlabels):
    label = equivalences.find(i)
    if renumber[label] is None:
      renumber[label] = next_label
      renumber[i] = next_label
      next_label += 1
    else:
      renumber[i] = renumber[label]

  if len(renumber):
    renumber = np.array(renumber)
    out[:] = renumber[out.reshape((sx*sy,))].reshape((sx,sy), order='C')

  return out, next_label

