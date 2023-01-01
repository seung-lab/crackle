import crackle
import numpy as np

import compresso

import pytest

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

def test_compress_decompress_z_range():
  labels = compresso.load("connectomics.npy.cpso.gz")

  x,y,z = tuple(np.random.randint(128,384, size=(3,)))
  print(x,y,z)
  cutout = labels[x:x+128,y:y+128,z:z+16]
  
  binary = crackle.compress(cutout)
  arr = crackle.CrackleArray(binary)

  recovered = arr[2:100,5:83,5:7]
  assert np.all(cutout[2:100,5:83,5:7] == recovered)

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

