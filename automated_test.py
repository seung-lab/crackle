import crackle
import numpy as np

import compresso

import pytest

DTYPE = [
  np.uint8, np.uint16, np.uint32, np.uint64,
  # np.int8,  np.int16,  np.int32,  np.int64
]

@pytest.mark.parametrize('dtype', DTYPE)
@pytest.mark.parametrize('markov_model_order', [0,1,2,3])
def test_compress_decompress_random(dtype, markov_model_order):
  labels = np.random.randint(0,5,size=(4,4,1), dtype=dtype)
  binary = crackle.compress(labels, markov_model_order=markov_model_order)
  recovered = crackle.decompress(binary)
  print(labels.T)
  print(recovered.T)
  assert np.all(labels == recovered)

  labels = np.random.randint(0,40,size=(1000,999,1), dtype=dtype)
  binary = crackle.compress(labels, markov_model_order=markov_model_order)
  recovered = crackle.decompress(binary)
  assert np.all(labels == recovered)

  labels = np.random.randint(0,40,size=(100,100,10), dtype=dtype)
  binary = crackle.compress(labels, allow_pins=True, markov_model_order=markov_model_order)
  recovered = crackle.decompress(binary)
  assert np.all(labels == recovered)

  labels = np.random.randint(0,40,size=(100,100,100), dtype=dtype)
  binary = crackle.compress(labels, markov_model_order=markov_model_order)
  recovered = crackle.decompress(binary)
  assert np.all(labels == recovered)

@pytest.mark.parametrize('allow_pins', [False,True])
@pytest.mark.parametrize('dtype', DTYPE)
def test_compress_decompress(dtype, allow_pins):
  labels = compresso.load("connectomics.npy.cpso.gz")

  x,y,z = tuple(np.random.randint(128,384, size=(3,)))
  print(x,y,z)
  cutout = labels[x:x+128,y:y+128,z:z+16].astype(dtype)
  
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

def test_connectomics_npy():
  labels = compresso.load("connectomics.npy.cpso.gz")
  binary = crackle.compress(labels, allow_pins=False)
  recovered = crackle.decompress(binary)
  assert np.all(labels == recovered)

def test_watershed():
  labels = compresso.load("ws.npy.cpso.gz")
  binary = crackle.compress(labels, allow_pins=False)
  recovered = crackle.decompress(binary)
  assert np.all(labels == recovered)

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

def test_num_labels():
  labels = np.random.randint(0,100, size=(100,100,10), dtype=np.uint32)
  binary = crackle.compress(labels, allow_pins=False)
  num_labels_1 = np.unique(labels).size  
  num_labels_2 = crackle.num_labels(binary)
  assert num_labels_1 == num_labels_2

  labels = np.ones((20,20,20), dtype=np.uint32)
  labels[2,2,2:5] = 4
  binary = crackle.compress(labels, allow_pins=True)
  num_labels_1 = np.unique(labels).size  
  num_labels_2 = crackle.num_labels(binary)
  assert num_labels_1 == num_labels_2

def test_contains():
  labels = np.random.randint(0,4, size=(100,100,10), dtype=np.uint32)
  binary = crackle.compress(labels)

  assert crackle.contains(binary, -1) == False
  assert crackle.contains(binary, 0) == True
  assert crackle.contains(binary, 1) == True
  assert crackle.contains(binary, 2) == True
  assert crackle.contains(binary, 3) == True
  assert crackle.contains(binary, 4) == False
  assert crackle.contains(binary, 21954124) == False
  assert crackle.contains(binary, -12998124) == False

def test_remap():
  labels = np.arange(0,1000).reshape((10,10,10)).astype(np.uint32)
  binary = crackle.compress(labels)

  labels = crackle.labels(binary)
  assert np.all(labels == np.arange(0,1000))

  binary = crackle.remap(binary, { i:i+1 for i in range(1000) }, preserve_missing_labels=False)
  labels = crackle.labels(binary)
  assert np.all(labels == np.arange(1,1001))

