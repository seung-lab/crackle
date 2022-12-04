# Crackle: Next gen. 3D segmentation compression codec.

```python
import crackle
import numpy

labels = np.load(...) # a 2D or 3D dense segmentation

binary = crackle.compress(labels)
labels = crackle.decompress(binary)

# get unique labels without decompressing
uniq = crackle.labels(binary) 

# Remap labels without decompressing. Could
# be useful for e.g. proofreading.
remapped = crackle.remap(
  binary, { 1: 2, 2: 3, ... },
  preserve_missing_labels=True
)

# for working with files
# if .gz is appended to the filename, the file will be
# automatically gzipped (or ungzipped)
crackle.save(labels, "example.ckl.gz")
labels = crackle.load("example.ckl.gz")
```

*This repository is currently highly experimental.*

Crackle is a new codec inspired by Compresso \[1\] for creating highly compressed 3D dense segmentation images. Compresso innovated by separating labels from boundary structures. There were conceptually four (but really five) elements in the format: header, labels, bit packed and RLE encoded binary image boundaries, and indeterminate boundary locations. 

Crackle improves upon Compresso by replacing the bit-packed boundary map with a "crack code" and also uses 3D information to reduce redundancy in labels using "pins". Like Compresso, Crackle uses a two pass compression strategy where the output of crackle may be further comrpessed with a bitstream compressor like gzip, bzip2, zstd, or lzma.

Based on preliminary experiments, it seems likely that the output of Crackle will be in the ballpark of 2x to 4x smaller than Compresso. The second stage compressed Crackle file will likely be about 70% the size of the equivalent Compresso file.

## Boundary Structure: Crack Code

Our different approach is partially inspired by the work of Zingaretti et al. \[2\]. We represent the boundary not by border voxels, but by a "crack code" that represents the edges between voxels. This code can be thought of as directions to draw edges on a graph where the vertices are where the corners of four pixels touch and the edges are the cracks in between them. 

Since this regular graph is 4-connected, each "move" in a cardinal direction can be described using two bits. To represent special symbols such as "branch" and "terminate", an impossible set of instructions on an undirected graph such as "left-right" or "up-down" can be used (occupying 4 bits). In order to avoid creating palendromic sequences such as (3, 0, 3) meaning (down, branch) but can be read (terminate, down), we can use the left-right impossible directions to rewrite it as (3, 2, 1).

While the image is 3D, we treat the image in layers because working in 3D introduces a large increase in geometric complexity (a cube has 6 faces, 12 edges, and 8 corners while a square has 4 edges and 4 corners). This increase in complexity would inflate the size of the crack code and make the implementation more difficult.

## Label Map: Method of Pins

Each 2D CCL region must has a label assigned. Due to the 2D nature of the crack code, we cannot use 3D CCL. However, for example, a solid cube of height 100 would need 100 labels to represent the same color on every slice as in Compresso.

It is still possible to reduce the amount of redundant information even without 3D CCL. For each label, we find a set of vertical line segments ("pins") that fully cover the label's 2D CCL regions. Sharp readers may note that this is the NP-hard set cover problem.

Once a reasonably small or minimal set of pins are found, they can be encoded in two forms:

Condensed Form: `[label][num_pins][pin_1][pin_2]...[pin_N]`
Fixed Width Form: `[label][pin_1][label][pin_2]...[label][pin_N]`
Pin Format: `[linear index of pin top][number of voxels to bottom]`

Fixed width example with label 1 with a pin between (1,1,1) and (1,1,5) on a 10x10x10 image: `[1][111][4]`

An alternative formulation `[label][idx1][idx2]` was shown in an experiment on `connectomics.npy.cpso` to compress slightly worse than Compresso labels. However, this alternative formulation theoretically allows arbitrary pin orientations and so might be useful for reducing the overall number of pins.

The condensed format is a bit smaller than the fixed width format, but the fixed width format enables rapid searches if the set of pins are sorted by either the label (enables fast `label in file`) or the likely more useful sorting by top index to filter candidate pins when performing random access to a z-slice.

## References


1. Matejek, B., Haehn, D., Lekschas, F., Mitzenmacher, M., Pfister, H., 2017. Compresso: Efficient Compression of Segmentation Data for Connectomics, in: Descoteaux, M., Maier-Hein, L., Franz, A., Jannin, P., Collins, D.L., Duchesne, S. (Eds.), Medical Image Computing and Computer Assisted Intervention − MICCAI 2017, Lecture Notes in Computer Science. Springer International Publishing, Cham, pp. 781–788. https://doi.org/10.1007/978-3-319-66182-7_89

2. Zingaretti, P., Gasparroni, M., Vecci, L., 1998. Fast chain coding of region boundaries. IEEE Transactions on Pattern Analysis and Machine Intelligence 20, 407–415. https://doi.org/10.1109/34.677272




