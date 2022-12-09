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

  codes = crackle.crackcode.create_crack_codes(labels, permissible=False)
  ans = [[9, 1, 0, 2, 3, 3, 0]]
  assert codes == ans

@pytest.mark.parametrize("size", [16, 64, 128, 1024])
def test_packed_encoding(size):
  chains = [[9, 0, 3, 1, 0, 2, 3, 3, 0, 1, 2]]
  ans = { 9: ['b', 1, 0, 2, 3, 't', 't'] }

  sx = sy = size

  packed_code = crackle.crackcode.pack_codes(chains, sx, sy)
  recovered = crackle.crackcode.unpack_binary(packed_code, sx, sy)

  assert ans == recovered

def test_remove_initial_branch():
  code = []
  assert crackle.crackcode.remove_initial_branch(code, 5, 5) == []

  code = [10, 'd','u','t']
  assert crackle.crackcode.remove_initial_branch(code, 5, 5) == code

  code = [10, 'b', 'd', 'r', 't']
  revised = crackle.crackcode.remove_initial_branch(code, 5, 5)
  assert revised == [17, 's', 'u', 'l', 's']

  code = [10, 'b', 'd', 'b', 'r', 't']
  revised = crackle.crackcode.remove_initial_branch(code, 5, 5)
  assert revised == [10, 'b', 'd', 'b', 'r', 't']


def test_compress_decompress_empty():
  labels = np.zeros((100,100,32), dtype=np.uint32)
  binary = crackle.compress(labels)
  recovered = crackle.decompress(binary)
  assert np.all(labels == recovered)

def test_compress_decompress_random():
  labels = np.random.randint(0,5,size=(4,4,1), dtype=np.uint32)
  binary = crackle.compress(labels)
  recovered = crackle.decompress(binary)
  assert np.all(labels == recovered)

  labels = np.random.randint(0,40,size=(1000,999,1), dtype=np.uint32)
  binary = crackle.compress(labels)
  recovered = crackle.decompress(binary)
  assert np.all(labels == recovered)

  labels = np.random.randint(0,40,size=(100,100,100), dtype=np.uint32)
  binary = crackle.compress(labels)
  recovered = crackle.decompress(binary)
  assert np.all(labels == recovered)

@pytest.mark.parametrize('i', range(3))
def test_compress_decompress(i):
  labels = compresso.load("connectomics.npy.cpso.gz")

  x,y,z = tuple(np.random.randint(128,384, size=(3,)))
  print(x,y,z)
  cutout = labels[x:x+128,y:y+128,z:z+16]
  
  binary = crackle.compress(cutout)
  recovered = crackle.decompress(binary)

  assert np.all(cutout == recovered)

def test_labels():
  labels = np.random.randint(0,100, size=(100,100,10), dtype=np.uint32)
  binary = crackle.compress(labels)
  uniq = np.unique(labels)  
  uniq2 = crackle.labels(binary)

  assert np.all(uniq == uniq2)

  labels = compresso.load("connectomics.npy.cpso.gz")
  x,y,z = tuple(np.random.randint(128,384, size=(3,)))
  print(x,y,z)
  labels = labels[x:x+128,y:y+128,z:z+16]
  binary = crackle.compress(labels)
  uniq = np.unique(labels)  
  uniq2 = crackle.labels(binary)
  assert np.all(uniq == uniq2)

def test_remap():
  labels = np.arange(0,1000).reshape((10,10,10))
  binary = crackle.compress(labels)

  labels = crackle.labels(binary)
  assert np.all(labels == np.arange(0,1000))

  binary = crackle.remap(binary, { i:i+1 for i in range(1000) }, preserve_missing_labels=False)
  labels = crackle.labels(binary)
  assert np.all(labels == np.arange(1,1001))  









