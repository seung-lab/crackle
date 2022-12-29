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
    
    s = time.time()
    ckl_binary = crackle.compress(cutout, force_flat=True)
    compress_time = time.time() - s

    s = time.time()
    cpso_binary = compresso.compress(cutout)
    cpso_time = time.time() - s

    raw_binary = cutout.tobytes('F')

    s = time.time()
    cckl_binary = zlib.compress(ckl_binary)
    compress_time_gz = compress_time + (time.time() - s)

    ccpso_binary = zlib.compress(cpso_binary)
    craw_binary = zlib.compress(raw_binary)

    mvxs = cutout.size / compress_time / 1e6
    mvxs_gz = cutout.size / compress_time_gz / 1e6
    mvxs_cpso = cutout.size / cpso_time / 1e6

    print(f"""
      ckl:     {len(ckl_binary)}    ({mvxs:.2f} MVx/sec)
      cpso:    {len(cpso_binary)}   ({len(ckl_binary)/len(cpso_binary)*100:.2f}%) ({mvxs_cpso:.2f} MVx/sec)
      raw:     {len(raw_binary)}   ({len(ckl_binary)/len(raw_binary)*100:.2f}%)
      ckl.gz:  {len(cckl_binary)}   ({mvxs_gz:.2f} MVx/sec)
      cpso.gz  {len(ccpso_binary)}  ({len(cckl_binary)/len(ccpso_binary)*100:.2f}%)
      raw:.gz  {len(craw_binary)}  ({len(cckl_binary)/len(craw_binary)*100:.2f}%)
    """, flush=True)

N = 3
shape = (128,128,64)

print("shape:", shape)

print("PINKY40 CUTOUTS (connectomics.npy)")
labels = compresso.load("benchmarks/connectomics.npy.cpso.gz")
run_sample(labels, shape, N)

print("WATERSHED CUTOUTS (ws.npy)")
labels = compresso.load("benchmarks/ws.npy.cpso.gz")
run_sample(labels, shape, N)

print("RANDOM NOISE [0,2000) uint32")
labels = np.random.randint(0,2000, size=(512,512,512), dtype=np.uint32)
run_sample(labels, shape, N)

print("EMPTY")
labels = np.zeros((512,512,512), dtype=np.uint32)
run_sample(labels, shape, 1)

print("SOLID ONES")
labels = np.ones((512,512,512), dtype=np.uint32)
run_sample(labels, shape, 1)

