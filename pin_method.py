from collections import defaultdict

import compresso
import numpy as np
import cc3d
from tqdm import tqdm

labels = compresso.load("connectomics.npy.cpso")
sx,sy,sz = labels.shape

cc_labels = np.zeros((sx,sy,sz), dtype=np.uint32)
for z in range(sz):
	cc_labels[:,:,z] = cc3d.connected_components(labels[:,:,z], connectivity=4)

pinsets = defaultdict(list)

for x in tqdm(range(sx)):
	for y in range(sy):
		label = labels[x,y,0]
		label_set = set([ cc_labels[x,y,0] ])
		z_start = 0
		for z in range(1, sz):
			cur = labels[x,y,z]
			if label != cur:
				pinsets[label].append(
					((x,y,z_start), (x,y,z), label_set)
				)
				label = cur
				z_start = z
				label_set = set()
			label_set |= set([ cc_labels[x,y,z] ])
		pinsets[label].append(
			((x,y,z_start), (x,y,z), label_set)
		)

def find_optimal_pins(pinset):
	sets = [ x[2] for x in pinset ]
	isets = [ [i,x] for i,x in enumerate(sets) ]

	universe = set.union(*sets)

	final_pins = []
	while len(universe):
		sizes = [ len(x[1]) for x in isets ]
		idx = sizes.index(max(sizes))
		i, cur = isets.pop(idx)
		universe -= cur
		for j, otherset in isets:
			otherset -= cur
		final_pins.append(pinset[i][:2])

	return final_pins

all_pins = {}

for label in tqdm(pinsets):
	pin = find_optimal_pins(pinsets[label])
	all_pins[label] = pin

npins = sum([ len(x) for x in all_pins.values() ])
nlabels = len(all_pins)

cpso = compresso.compress(labels)
raws = compresso.raw_labels(cpso).tobytes()

toidx = lambda tup: int(tup[0] + sx * tup[1] + sx * sy * tup[2])

ZBYTES = 2 if sz > 127 else 1

def condensed_binary(all_pins):
	linear = []
	for label, pins in all_pins.items():
		# linear.append(int(len(pins)).to_bytes(1, 'little'))
		for pin in pins:
			linear.append(int(label).to_bytes(4, 'little'))
			linear.append(toidx(pin[0]).to_bytes(4, 'little'))
			linear.append(int(pin[1][2] - pin[0][2]).to_bytes(ZBYTES, 'little'))
	return b''.join(linear)

def verbose_binary(all_pins):
	linear = []
	for label, pins in all_pins.items():
		linear.append(int(label).to_bytes(4, 'little'))
		linear.append(int(len(pins)).to_bytes(1, 'little'))
		for pin in pins:
			linear.append(toidx(pin[0]).to_bytes(4, 'little'))
			linear.append(int(pin[1][2] - pin[0][2]).to_bytes(ZBYTES, 'little'))
	return b''.join(linear)


# pcost = 2 * 4 * npins + (4+1) * nlabels
# pcost = 3 * 4 * npins

pinbin = condensed_binary(all_pins)
pcost = len(pinbin)
ccost = len(raws)

print(f"shape: {labels.shape}")
print(f"pin cost: {pcost} bytes")
print(f"original cost: {ccost} bytes")

print(f"Relative Size (pins/orig): {pcost/ccost*100:0.1f}%")

with open("pins.raw", "wb") as f:
	f.write(pinbin)

with open("orig.raw", "wb") as f:
	f.write(raws)












