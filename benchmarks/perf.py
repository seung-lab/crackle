import crackle
import compresso

import numpy as np
import deflate

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
    gzip_bin = deflate.gzip_compress(ckl_binary, 6)
    gzip_time = time.time() - s

    s = time.time()
    deflate.gzip_decompress(gzip_bin)
    gunzip_time = time.time() - s

    s = time.time()
    crackle.decompress(ckl_binary)
    decompress_time = time.time() - s    

    mvxs = lambda t: cutout.size / t / 1e6

    print(f"""
      compress     :  {mvxs(compress_time):.2f} MVx/sec ({len(ckl_binary)} bytes, {len(ckl_binary)/cutout.nbytes*100:.1f}%)
      compress+gz  :  {mvxs(compress_time+gzip_time):.2f} MVx/sec ({len(gzip_bin)} bytes, {len(gzip_bin)/len(ckl_binary)*100:.1f}%)
      decompress   :  {mvxs(decompress_time):.2f} MVx/sec
      decompress+gz:  {mvxs(decompress_time+gunzip_time):.2f} MVx/sec
    """, flush=True)

N = 3
shape = (256,256,64)

print("shape:", shape)

print("PINKY40 CUTOUTS (connectomics.npy)")
labels = compresso.load("benchmarks/connectomics.npy.cpso.gz")
run_sample(labels, shape, N)

print("WATERSHED CUTOUTS (ws.npy)")
labels = compresso.load("benchmarks/ws.npy.cpso.gz")
run_sample(labels, shape, N)

print("RANDOM NOISE [0,2000) uint32")
labels = np.random.randint(0,2000, size=(512,512,512), dtype=np.uint32)
labels = np.asfortranarray(labels)
run_sample(labels, shape, N)

print("BINARY NOISE [0,1] uint8 (pathological case)")
labels = np.random.randint(0,2, size=(512,512,512), dtype=np.uint8)
labels = np.asfortranarray(labels)
run_sample(labels, shape, N)

print("EMPTY")
labels = np.zeros((512,512,512), dtype=np.uint32, order="F")
run_sample(labels, shape, 1)

print("SOLID ONES")
labels = np.ones((512,512,512), dtype=np.uint32, order="F")
run_sample(labels, shape, 1)

