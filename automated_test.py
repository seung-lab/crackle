import crackle
import crackle.crackcode
import numpy as np

import compresso

import pytest

def test_create_graph():
  labels = np.zeros((3,3), dtype=np.uint32)
  labels[1,1] = 1

  G = crackle.crackcode.create_graph(labels)

  assert set(G.nodes) == set([5,9,6,10])
  assert set(G.edges) == set([(5,9), (5,6), (6,10), (9,10)])

def test_create_crack_codes():
  labels = np.zeros((3,3), dtype=np.uint32)
  labels[1,1] = 1

  codes = crackle.crackcode.create_crack_codes(labels)
  ans = [[9, 0, 3, 1, 0, 2, 3, 3, 0, 1, 2]]
  assert codes == ans

def test_packed_encoding():
  chains = [[9, 0, 3, 1, 0, 2, 3, 3, 0, 1, 2]]
  ans = { 9: ['b', 1, 0, 2, 3, 't', 't'] }
  packed_code = crackle.crackcode.pack_codes(chains)
  recovered = crackle.crackcode.unpack_binary(packed_code)

  assert ans == recovered

def test_compress_decompress_empty():
  labels = np.zeros((100,100,32), dtype=np.uint32)
  binary = crackle.compress(labels)
  recovered = crackle.decompress(binary)
  assert np.all(labels == recovered)

@pytest.mark.parametrize('i', range(10))
def test_compress_decompress(i):
  labels = compresso.load("connectomics.npy.cpso.gz")

  x,y,z = tuple(np.random.randint(128,384, size=(3,)))
  print(x,y,z)
  cutout = labels[x:x+128,y:y+128,z:z+16]
  
  binary = crackle.compress(cutout)
  recovered = crackle.decompress(binary)

  assert np.all(cutout == recovered)





