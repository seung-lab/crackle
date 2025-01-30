import random

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

  labels = np.random.randint(0,40,size=(256,255,1), dtype=dtype)
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

@pytest.mark.parametrize('label', [9999999999, 63408621, 63408621, 28792756])
def test_decompress_binary_label_flat(label):
  labels = compresso.load("connectomics.npy.cpso.gz")
  binary = crackle.compress(labels)

  recovered = crackle.decompress(binary)
  binimg1 = recovered == label
  binimg2 = crackle.decompress(binary, label)

  assert np.all(labels == recovered)
  assert np.all(binimg1 == binimg2)

@pytest.mark.parametrize('allow_pins', [True])
def test_decompress_binary_label_random(allow_pins):
  labels = compresso.load("connectomics.npy.cpso.gz")

  x,y,z = tuple(np.random.randint(128,256, size=(3,)))
  print(x,y,z)
  labels = np.asfortranarray(labels[x:x+256,y:y+256,z:z+32])

  binary = crackle.compress(labels, allow_pins=allow_pins)

  label = random.choice(np.unique(labels))

  recovered = crackle.decompress(binary)
  binimg1 = recovered == label
  binimg2 = crackle.decompress(binary, label)

  assert np.all(labels == recovered)
  assert np.all(binimg1 == binimg2)

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

  shape = [100,100,50]
  labels = np.arange(100*100*50, dtype=np.uint32).reshape(shape, order="F")
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

  crackle.decompress(binary) # validate crcs

  binary = crackle.remap(binary, { i:i+1 for i in range(1002) }, preserve_missing_labels=True)
  labels = crackle.labels(binary)
  assert np.all(labels == np.arange(2,1002))

  ans = np.arange(2,1002)
  ans[:48] += 1

  binary = crackle.remap(binary, { i:i+1 for i in range(50) }, preserve_missing_labels=True)
  labels = crackle.labels(binary)
  assert np.all(labels == ans)

def test_big_remap():
  import time
  labels = np.arange(0,10000000).reshape((100,100,1000)).astype(np.uint32)
  binary = crackle.compress(labels)
  mapping = { i:i+1 for i in range(10000000) }
  s = time.time()
  crackle.remap(binary, mapping, preserve_missing_labels=False, in_place=True, parallel=0)  
  print(f"{time.time() - s:.3f}")

  labels = crackle.labels(binary)
  assert np.all(labels == np.arange(1,10000001))

def test_remap_sorted():
  labels = np.random.randint(0,7, size=(10,10,10)).astype(np.uint32)
  binary = crackle.compress(labels)

  remap = {0:10, 1:40, 2:20, 3:90, 4:25, 5:70, 6:90}

  binary = crackle.remap(binary, remap, preserve_missing_labels=True)
  new_labels = list(remap.values())
  new_labels.sort()
  remapped_binary_labels = list(crackle.labels(binary))
  print(remapped_binary_labels)
  remapped_binary_labels.sort()

  assert remapped_binary_labels == new_labels

  for label in new_labels:
    assert crackle.contains(binary, label)

@pytest.mark.parametrize("allow_pins", [False, True])
def test_zstack_ones(allow_pins):

  sz = 100

  binary = crackle.zstack([
    np.ones([512,512], dtype=np.uint32)
    for z in range(sz)
  ])

  head = crackle.header(binary)
  assert head.sx == 512
  assert head.sy == 512
  assert head.sz == sz

  assert head.data_width == 4
  assert head.stored_data_width == 1

  assert crackle.labels(binary) == [1]

  labels = np.ones([512,512,sz], dtype=np.uint32)
  binary2 = crackle.compress(labels, allow_pins=allow_pins)

  if not allow_pins:
    assert binary == binary2

  arr = crackle.decompress(binary)
  assert np.all(arr == 1)

  sz = 100

  binary = crackle.zstack([
    np.ones([512,512,2], dtype=np.uint32)
    for z in range(sz // 2)
  ])
  labels = np.ones([512,512,sz], dtype=np.uint32)
  binary2 = crackle.compress(labels, allow_pins=allow_pins)

  if not allow_pins:
    assert binary == binary2

  sz = 10

  image = np.random.randint(0,255, size=(64, 64, sz), dtype=np.uint8)
  binary = crackle.zstack([
    image[:,:,2*z:2*z+2]
    for z in range(sz // 2)
  ])

  recovered = crackle.decompress(binary)
  assert np.all(image == recovered)

@pytest.mark.parametrize('allow_pins', [False,True])  
def test_zstack_real_data(allow_pins):
  labels = compresso.load("connectomics.npy.cpso.gz")
  binaries = []
  for z in range(0, 512 // 64):
    slab = labels[:,:,z*64:(z+1)*64]
    binary = crackle.compress(slab, allow_pins=allow_pins)
    binaries.append(binary)

  binary = crackle.zstack(binaries)
  recovered = crackle.decompress(binary)

  assert np.all(labels == recovered)

def test_zsplit():
  labels = np.ones([256,256,7], dtype=np.uint32)
  binary = crackle.compress(labels)
  (before, middle, after) = crackle.zsplit(binary, 3)

  before = crackle.decompress(before)
  assert np.all(before == labels[:,:,:3])

  middle = crackle.decompress(middle)
  assert np.all(middle == labels[:,:,3:4])

  after = crackle.decompress(after)
  assert np.all(after == labels[:,:,4:])

def test_zshatter():
  labels = np.ones([256,256,7], dtype=np.uint32)
  binary = crackle.compress(labels)
  levels = crackle.zshatter(binary)

  for i in range(7):
    arr = crackle.decompress(levels[i])
    assert np.all(arr == 1)


def test_edit_array():
  labels = np.ones([256,256,10], dtype=np.uint32)
  binary = crackle.compress(labels)
  arr = crackle.CrackleArray(binary)

  res = np.full(labels.shape, 2, dtype=np.uint32)
  arr[:,:,:] = res
  assert np.all(arr[:] == res)

  labels = np.ones([256,256,10], dtype=np.uint32)
  binary = crackle.compress(labels)
  arr = crackle.CrackleArray(binary)

  res = np.full([256, 256, 5], 2, dtype=np.uint32)
  arr[:,:,3:8] = res

  ans = np.ones(labels.shape, dtype=np.uint32)
  ans[:,:,3:8] = 2

  assert np.all(arr[:] == ans)

def test_full():
  shape = [100,200,50]
  binary = crackle.full(shape, 0, dtype=np.uint16, order="F")
  ans = np.full(shape, 0, dtype=np.uint16, order="F")

  binary2 = crackle.compress(ans)

  assert binary == binary2

  arr = crackle.decompress(binary)

  assert np.all(arr == ans)

def test_zeros():
  binary = crackle.zeros([10,10,10], dtype=np.uint64, order="F")
  arr = crackle.decompress(binary)
  ans = np.zeros([10,10,10], dtype=np.uint64, order="F")
  assert crackle.labels(binary) == [ 0 ]
  assert np.all(arr == ans)

@pytest.mark.parametrize("scalar", [0, 1, 5, 255, 70000])
def test_add(scalar):
  binary = crackle.zeros([10,10,10], dtype=np.uint64, order="F")
  print(crackle.components(binary))

  arr = crackle.CrackleArray(binary)
  assert arr.labels() == [ 0 ]
  arr = arr + scalar
  assert arr.labels() == [ scalar ]

  assert np.all(
    arr.decompress() == np.full([10,10,10], scalar, dtype=np.uint64, order="F")
  )


def test_refit():
  binary = crackle.zeros([10,10,10], dtype=np.uint64, order="F")
  arr = crackle.CrackleArray(binary)
  arr = arr.refit()
  assert arr.dtype == np.uint8

def test_contiguous_fortran():
  binary = crackle.zeros([10,10,10], dtype=np.uint64, order="F")
  
  binary2 = crackle.asfortranarray(binary)
  assert binary == binary2

  binary2 = crackle.ascontiguousarray(binary)
  assert binary != binary2

  head = crackle.header(binary2)
  assert head.fortran_order == False

  arr = crackle.decompress(binary2)
  assert arr.flags.f_contiguous == False
  assert arr.flags.c_contiguous == True

  binary2 = crackle.asfortranarray(binary2)
  assert binary2 == binary

def test_point_cloud():
  input_arr = np.array([
    [0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0],
  ], dtype=np.uint32, order="F")

  binary = crackle.compress(input_arr)
  ptc = crackle.point_cloud(binary, 0)[:,:2]

  sx, sy = input_arr.shape

  result = []
  for x in range(sx):
    result.append([x,0])
    result.append([x,sy-1])
  for y in range(1, sy -1):
    result.append([0,y])
    result.append([sx-1,y])

  result.append([0,0]) # extra copy of 0,0

  result = np.array(result)
  result = np.sort(result, axis=0)
  ptc = np.sort(ptc, axis=0)

  assert np.all(result == ptc)

@pytest.mark.parametrize("allow_pins", [False, True])
@pytest.mark.parametrize("dtype", [np.uint16, np.uint32])
def test_min_max(allow_pins, dtype):
  labels = np.random.randint(0,500,size=(20,20,20), dtype=dtype)
  binary = crackle.compress(labels, allow_pins=allow_pins)

  ans_min = np.min(labels)
  ans_max = np.max(labels)

  assert ans_min == crackle.min(binary)
  assert ans_max == crackle.max(binary)

  labels = np.ones((20,20,20), dtype=dtype)
  binary = crackle.compress(labels, allow_pins=allow_pins)

  ans_min = np.min(labels)
  ans_max = np.max(labels)

  assert ans_min == crackle.min(binary)
  assert ans_max == crackle.max(binary)

@pytest.mark.parametrize("dtype", [np.uint8, np.uint16, np.uint32, np.uint64])
def test_crc_check_one_bit_error(dtype):
  arr = np.random.randint(0, 255, size=[11,23,37], dtype=dtype)
  binary = bytearray(crackle.compress(arr))

  head = crackle.header(binary) # crashes if bad
  binary[0] = ord('z')

  try:
    head = crackle.header(binary) # crashes if bad
    assert False
  except:
    pass

  binary[0] = ord('c') # fixes
  head = crackle.header(binary) # crashes if bad
  binary[4] = 100
  
  try:
    head = crackle.header(binary) # crashes if bad
    assert False
  except:
    pass

  binary[4] = 1
  head = crackle.header(binary) # crashes if bad

  # Hamming Distance 1
  # now to see if crc actually works....
  # tests 28, 29 are false positives btw
  for i in range(5, 30):
    for j in range(8):
      val = binary[i]
      binary[i] |= (0b1 << j)

      try:
        head = crackle.header(binary) # crashes if bad
        assert False, i
      except:
        pass

      binary[i] = val

  for i in range(5, 30):
    for j in range(8):
      val = binary[i]
      binary[i] = binary[i] & ~(0b1 << j)

      try:
        head = crackle.header(binary) # crashes if bad
        assert False, i
      except:
        pass

      binary[i] = val


  # Hamming Distance 2
  orig_binary = bytearray(binary)
  N = 184

  for i in range(N):
    for j in range(i, N):
      binary = bytearray(orig_binary)

      ii = i // 8
      jj = j // 8
      iib = i - ii * 8
      jjb = j - jj * 8

      binary[ii] |= (0b1 << iib)
      binary[jj] |= (0b1 << jjb)

      try:
        head = crackle.header(binary) # crashes if bad
        assert False, i
      except:
        pass

  for i in range(N):
    for j in range(i, N):
      binary = bytearray(orig_binary)

      ii = i // 8
      jj = j // 8
      iib = i - ii * 8
      jjb = j - jj * 8

      binary[ii] = binary[ii] & ~(0b1 << iib)
      binary[jj] = binary[jj] & ~(0b1 << jjb)

      try:
        head = crackle.header(binary) # crashes if bad
        assert False, i
      except:
        pass

def test_array_works_z():
  labels = np.arange(0,64).reshape((4,4,4), order="F").astype(np.uint32)
  compressed = crackle.compress(labels)
  arr = crackle.CrackleArray(compressed)
  assert np.all(arr[:,:,3] == labels[:,:,3])

@pytest.mark.parametrize('allow_pins', [False,True])
def test_reencode(allow_pins):
  labels = compresso.load("connectomics.npy.cpso.gz")
  binary = crackle.compress(labels, allow_pins=allow_pins)
  
  markov_binary = crackle.codec.reencode(binary, markov_model_order=5)
  assert len(markov_binary) < len(binary)

  recovered = crackle.decompress(binary)
  assert np.all(recovered == labels)
  del recovered

  recovered = crackle.decompress(markov_binary)
  assert np.all(labels == recovered)  

def test_voxel_counts():
  labels = np.zeros([1,1,1])
  binary = crackle.compress(labels)
  vc_cts = crackle.voxel_counts(binary)
  assert vc_cts[0] == 1

  labels = compresso.load("connectomics.npy.cpso.gz")
  binary = crackle.compress(labels)

  vc_cts = crackle.voxel_counts(binary)

  uniq, cts = np.unique(labels, return_counts=True)
  cts_gt = { u:ct for u,ct in zip(uniq, cts) }
  assert vc_cts == cts_gt

  vc_ct = crackle.voxel_counts(binary, label=25024949)
  assert vc_ct == cts_gt[25024949]

  labels = np.zeros([100,101,103])
  binary = crackle.compress(labels)

  vc_cts = crackle.voxel_counts(binary)

  uniq, cts = np.unique(labels, return_counts=True)
  cts_gt = { u:ct for u,ct in zip(uniq, cts) }
  assert vc_cts == cts_gt


def test_spurious_branch_elimination():
  arr = np.array([
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,1,1,2,2,0,0,0,0],
    [0,0,1,1,2,2,0,0,0,0],
    [0,0,4,4,3,3,0,0,0,0],
    [0,0,4,4,3,3,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
  ], dtype=np.uint8).T

  binary = crackle.compress(arr)
  recovered = crackle.decompress(binary)[:,:,0]

  assert np.all(recovered == arr)

  arr = np.array([
    [  0, 139, 139, 139, 139],
    [  0, 139,   0, 139, 139],
    [  0, 161,   0,   0, 161],
    [161, 161, 161, 161, 161],
  ], dtype=np.uint8).T
  binary = crackle.compress(arr)
  recovered = crackle.decompress(binary)[:,:,0]

  assert np.all(recovered == arr)









