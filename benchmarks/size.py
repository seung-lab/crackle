import crackle
import compresso

import numpy as np
import zlib
import brotli

import time

def sample(labels, shape):
  x,y,z = tuple(np.random.randint(0,max(labels.shape) - max(shape), size=(3,)))
  return labels[x:x+shape[0],y:y+shape[1],z:z+shape[2]]

def run_sample(labels, shape, N):
  for i in range(N):
    cutout = sample(labels, shape)
    run(cutout)

def run(labels):
  ckl_binary = crackle.compress(labels, allow_pins=False)    
  mkv_binary = crackle.compress(labels, allow_pins=False, markov_model_order=5)
  pin_binary = crackle.compress(labels, allow_pins=True)
  cpso_binary = compresso.compress(labels)
  raw_binary = labels.tobytes('F')

  for method in ['','gz','br']:


    if method == '':
      fn = lambda x: x
      ext = ''
    elif method == 'gz':
      fn = zlib.compress
      ext = '.gz'
    elif method == 'br':
      fn = brotli.compress
      ext = '.br'

    cckl_binary = fn(ckl_binary)
    ccpso_binary = fn(cpso_binary)
    craw_binary = fn(raw_binary)
    cpin_binary = fn(pin_binary)
    cmkv_binary = fn(mkv_binary)

    print(f"""
      ckl{ext}:      {len(cckl_binary): 9}   ({len(cckl_binary)/len(craw_binary)*100:.2f}%)
      ckl{ext}(pin): {len(cpin_binary): 9}   ({len(cpin_binary)/len(craw_binary)*100:.2f}%) 
      ckl{ext}(mkv): {len(cmkv_binary): 9}   ({len(cmkv_binary)/len(craw_binary)*100:.2f}%)
      cpso{ext}:     {len(ccpso_binary): 9}  ({len(ccpso_binary)/len(craw_binary)*100:.2f}%)
      raw{ext}:      {len(craw_binary): 9}  ({len(craw_binary)/len(craw_binary)*100:.2f}%)""", flush=True)

N = 1
shape = (256,256,64)

print("shape:", shape)

print("PINKY40 CUTOUTS (connectomics.npy)")
labels = compresso.load("benchmarks/connectomics.npy.cpso.gz")
run_sample(labels, shape, N)

print("BINARY IMAGE (connectomics.npy, label 67699431)")
labels = compresso.load("benchmarks/connectomics.npy.cpso.gz")
labels = labels == 67699431
run(labels)

print("WATERSHED CUTOUTS (ws.npy)")
labels = compresso.load("benchmarks/ws.npy.cpso.gz")
run_sample(labels, shape, N)

print("EMPTY")
labels = np.zeros((512,512,512), dtype=np.uint32)
run_sample(labels, shape, 1)

print("SOLID ONES")
labels = np.ones((512,512,512), dtype=np.uint32)
run_sample(labels, shape, 1)

print("RANDOM NOISE [0,2000) uint32")
labels = np.random.randint(0,2000, size=(512,512,128), dtype=np.uint32)
run_sample(labels, shape, N)

print("BINARY NOISE [0,1] uint8 (pathological case)")
labels = np.random.randint(0,2, size=(512,512,128), dtype=np.uint8)
run_sample(labels, shape, N)
