import crackle
import compresso

import numpy as np
import zlib

import time

def sample(labels, shape):
  x,y,z = tuple(np.random.randint(0,max(labels.shape) - max(shape), size=(3,)))
  return labels[x:x+shape[0],y:y+shape[1],z:z+shape[2]]

def run_sample(labels, shape, N):
  for i in range(N):
    cutout = sample(labels, shape)
    run(cutout)

def run(labels):
    s = time.time()
    ckl_binary = crackle.compress(labels, allow_pins=False)
    compress_time = time.time() - s

    s = time.time()
    mkv_binary = crackle.compress(labels, allow_pins=False, markov_model_order=5)
    mkv_compress_time = time.time() - s

    s = time.time()
    pin_binary = crackle.compress(labels, allow_pins=True)
    pin_compress_time = time.time() - s

    s = time.time()
    cpso_binary = compresso.compress(labels)
    cpso_time = time.time() - s

    raw_binary = labels.tobytes('F')

    s = time.time()
    cckl_binary = zlib.compress(ckl_binary)
    compress_time_gz = compress_time + (time.time() - s)

    ccpso_binary = zlib.compress(cpso_binary)
    craw_binary = zlib.compress(raw_binary)
    cpin_binary = zlib.compress(pin_binary)
    cmkv_binary = zlib.compress(mkv_binary)

    print(f"""
      ckl:         {len(ckl_binary): 9}
      ckl(pin):    {len(pin_binary): 9}  ({len(pin_binary)/len(ckl_binary)*100:.2f}%) (vs flat)
      ckl(mkv):    {len(mkv_binary): 9}  ({len(mkv_binary)/len(ckl_binary)*100:.2f}%) (vs literals)
      cpso:        {len(cpso_binary): 9}  ({len(ckl_binary)/len(cpso_binary)*100:.2f}%) 
      raw:         {len(raw_binary): 9}  ({len(ckl_binary)/len(raw_binary)*100:.2f}%)
      ckl.gz:      {len(cckl_binary): 9}   
      ckl.gz(pin): {len(cpin_binary): 9}   ({len(cpin_binary)/len(cckl_binary)*100:.2f}%) (vs flat)
      ckl.gz(mkv): {len(cmkv_binary): 9}   ({len(cmkv_binary)/len(cckl_binary)*100:.2f}%) (vs literals)
      cpso.gz:     {len(ccpso_binary): 9}  ({len(cckl_binary)/len(ccpso_binary)*100:.2f}%)
      raw.gz:      {len(craw_binary): 9}  ({len(cckl_binary)/len(craw_binary)*100:.2f}%)
    """, flush=True)

N = 3
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

print("RANDOM NOISE [0,2000) uint32")
labels = np.random.randint(0,2000, size=(512,512,512), dtype=np.uint32)
run_sample(labels, shape, N)

print("BINARY NOISE [0,1] uint8 (pathological case)")
labels = np.random.randint(0,2, size=(512,512,512), dtype=np.uint8)
run_sample(labels, shape, N)

print("EMPTY")
labels = np.zeros((512,512,512), dtype=np.uint32)
run_sample(labels, shape, 1)

print("SOLID ONES")
labels = np.ones((512,512,512), dtype=np.uint32)
run_sample(labels, shape, 1)

