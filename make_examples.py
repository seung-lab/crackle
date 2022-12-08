import crackle
import compresso

import numpy as np
import gzip
import fastremap

labels = compresso.load("connectomics.npy.cpso.gz")

pts = [
[257, 352, 238],
[258, 135, 358],
[268, 366, 273],
[230, 244, 310],
[167, 213, 209],
[318, 265, 369],
[159, 187, 297],
[167, 263, 239],
[161, 361, 323],
[371, 298, 248]
]

# for i in range(10):
# 	x,y,z = tuple(np.random.randint(128,384, size=(3,)))
# 	print(x,y,z)

for i, (x,y,z) in enumerate(pts):
	cutout = labels[x:x+128,y:y+128,z:z+64]

	crackle.save(cutout, f"examples/test-{i}.ckl")
	compresso.save(cutout, f"examples/test-{i}.cpso")
