import crackle
import numpy as np

import compresso

import pytest

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

@pytest.mark.parametrize('allow_pins', [False,True])
@pytest.mark.parametrize('i', range(3))
def test_compress_decompress(i, allow_pins):
  labels = compresso.load("connectomics.npy.cpso.gz")

  x,y,z = tuple(np.random.randint(128,384, size=(3,)))
  print(x,y,z)
  cutout = labels[x:x+128,y:y+128,z:z+16]
  
  binary = crackle.compress(cutout, allow_pins=allow_pins)
  recovered = crackle.decompress(binary)

  assert np.all(cutout == recovered)

@pytest.mark.parametrize('allow_pins', [False,True])
def test_compress_decompress_z_range(allow_pins):
  labels = compresso.load("connectomics.npy.cpso.gz")

  x,y,z = tuple(np.random.randint(128,384, size=(3,)))
  print(x,y,z)
  cutout = labels[x:x+128,y:y+128,z:z+16]
  
  binary = crackle.compress(cutout, allow_pins=allow_pins)
  arr = crackle.CrackleArray(binary)

  recovered = arr[2:100,5:83,5:7]
  assert np.all(cutout[2:100,5:83,5:7] == recovered)

  recovered = arr[2:100,5:83,0]
  assert np.all(cutout[2:100,5:83,0] == recovered)

  recovered = arr[2:100,5:83,-1]
  assert np.all(cutout[2:100,5:83,-1] == recovered)

@pytest.mark.parametrize('allow_pins', [False,True])
def test_empty(allow_pins):
  labels = np.zeros((0,0,0), dtype=np.uint32, order="F")
  compressed = crackle.compress(labels, allow_pins=allow_pins)
  reconstituted = crackle.decompress(compressed)
  assert np.all(labels == reconstituted)
  assert np.all(np.unique(labels) == crackle.labels(compressed))

@pytest.mark.parametrize('allow_pins', [False,True])
def test_black(allow_pins):
  labels = np.zeros((100,100,100), dtype=np.uint32, order="F")
  compressed = crackle.compress(labels, allow_pins=allow_pins)
  reconstituted = crackle.decompress(compressed)
  assert np.all(labels == reconstituted)
  assert np.all(np.unique(labels) == crackle.labels(compressed))

@pytest.mark.parametrize('allow_pins', [False,True])
def test_uniform_field(allow_pins):
  labels = np.zeros((100,100,100), dtype=np.uint32, order="F") + 1
  compressed = crackle.compress(labels, allow_pins=allow_pins)
  reconstituted = crackle.decompress(compressed)
  assert len(compressed) < labels.nbytes
  assert np.all(labels == reconstituted)
  assert np.all(np.unique(labels) == crackle.labels(compressed))

  labels = np.zeros((100,100,100), dtype=np.uint32, order="F") + np.iinfo(np.uint32).max
  compressed2 = crackle.compress(labels, allow_pins=allow_pins)
  reconstituted = crackle.decompress(compressed2)
  assert len(compressed2) < labels.nbytes
  assert np.all(labels == reconstituted)
  assert np.all(np.unique(labels) == crackle.labels(compressed2))

@pytest.mark.parametrize('allow_pins', [False,True])
def test_arange_field(allow_pins):
  labels = np.arange(0,1024).reshape((16,16,4)).astype(np.uint32)
  compressed = crackle.compress(labels, allow_pins=allow_pins)
  reconstituted = crackle.decompress(compressed)
  assert np.all(labels == reconstituted)
  assert np.all(np.unique(labels) == crackle.labels(compressed))

  labels = np.arange(1,1025).reshape((16,16,4)).astype(np.uint32)
  compressed = crackle.compress(labels, allow_pins=allow_pins)
  reconstituted = crackle.decompress(compressed)
  assert np.all(labels == reconstituted)
  assert np.all(np.unique(labels) == crackle.labels(compressed))

@pytest.mark.parametrize('allow_pins', [False,True])
def test_2d_arange_field(allow_pins):
  labels = np.arange(0,16*16).reshape((16,16,1)).astype(np.uint32)
  compressed = crackle.compress(labels, allow_pins=allow_pins)
  reconstituted = crackle.decompress(compressed)
  assert np.all(labels == reconstituted)
  assert np.all(np.unique(labels) == crackle.labels(compressed))

@pytest.mark.parametrize('allow_pins', [False,True])
def test_2_field(allow_pins):
  labels = np.arange(0,1024).reshape((16,16,4)).astype(np.uint32)
  compressed = crackle.compress(labels, allow_pins=allow_pins)
  reconstituted = crackle.decompress(compressed)
  assert np.all(labels == reconstituted)
  assert np.all(np.unique(labels) == crackle.labels(compressed))
  
  labels[2,2,1] = np.iinfo(np.uint32).max
  compressed = crackle.compress(labels, allow_pins=allow_pins)
  reconstituted = crackle.decompress(compressed)
  assert np.all(labels == reconstituted)
  assert np.all(np.unique(labels) == crackle.labels(compressed))

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

