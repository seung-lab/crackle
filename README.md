# Crackle: Next gen. 3D segmentation compression codec.

```bash
# Command Line Interface
crackle data.npy # creates data.ckl
crackle -d data.ckl # recovers data.npy
```

```python
import crackle
import numpy

labels = np.load("example.npy") # a 2D or 3D dense segmentation

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

arr = crackle.CrackleArray(binary)
res = arr[:10,:10,:10] # array slicing (efficient z ranges)
20 in arr # highly efficient search
```

*This repository is currently experimental.*

Crackle is a new codec inspired by Compresso \[1\] for creating highly compressed 3D dense segmentation images. Compresso innovated by separating labels from boundary structures. There were conceptually four (but really five) elements in the format: header, labels, bit packed and RLE encoded binary image boundaries, and indeterminate boundary locations. 

Crackle improves upon Compresso by replacing the bit-packed boundary map with a "crack code" and also uses 3D information to reduce redundancy in labels using "pins". Like Compresso, Crackle uses a two pass compression strategy where the output of crackle may be further comrpessed with a bitstream compressor like gzip, bzip2, zstd, or lzma.

Based on benchmarks, it seems likely that the output of Crackle will be in the ballpark of 20\% to 50\% the size of Compresso. The second stage compressed Crackle file will likely be about 60\% to 85\% the size of the equivalent Compresso file.


## Versions

| Major Version | Format Version | Description                                                    |
|---------------|----------------|----------------------------------------------------------------|
| 1             | 0              | Still experimental. |


## Stream Format

| Section   | Bytes                                     | Description                                                                                                     |
|-----------|-------------------------------------------|-----------------------------------------------------------------------------------------------------------------|
| Header    | 24                                        | Metadata incl. length of fields.                                                                                |
| Crack Index     | header.sz * sizeof(uint32)             | Number of bytes for the crack codes in each slice.
| Labels       | header.num_label_bytes        | Can be either "flat" labels or "pins". Describes how to color connected components.                                                                                   |
| Crack Codes    | Variable length.           | Instructions for drawing crack boundaries.             |

### Header


| Attribute         | Value             | Type    | Description                                     |
|-------------------|-------------------|---------|-------------------------------------------------|
| magic             | crkl              | char[4] | File magic number.                              |
| format_version    | 0                 | u8      | Stream version.                   |
| format_field      | bitfield          | u16     | See below.                 |
| sx, sy, sz        | >= 0              | u32 x 3 | Size of array dimensions.                       |
| grid_size         | log2(grid_size)   | u8      | Stores log2 of grid dimensions in voxels.          |
| num_label_bytes   | Any.              | u8      | Number of bytes of the labels section. Note the labels come in at least two format types.          |


Format Field (u8): DDSSCLLFGRRRRRRR (each letter represents a bit, left is LSB)

DD: 2^(DD) = byte width of returned array (1,2,4,8 bytes)  
SS: 2^(SS) = byte width of stored labels (sometimes you can store values in 2 bytes when the final array is 8 bytes)  
C: 1: crack codes denote impermissible boundaries 0: they denote permissible boundaries.  
LL: 0: "flat" label format, 1: fixed width pins (unused?) 2: variable width pins 3: reserved  
F: whether the array is to be rendered as C (0) or F (1) order
G: Signed (if (1), data are signed int, otherwise unsigned int)
R: Reserved

### Flat Label Format

`| NUM_LABELS (u64) | UNIQUE LABELS (NUM_LABELS \* STORED_DATA_WIDTH) | NUM CONNECTED COMPONENTS PER A GRID (sizeof(sx \* sy) \* sz) | INDEX INTO UNIQUE LABELS FOR EACH CCL ID (sizeof(NUM_LABELS) \* ... to end of section) |`

Flat labels are random access read, allow efficient reading of unique labels, efficient remapping, and efficient search for a given label's existence. Since the connected component labels can often use a smaller byte width than the unique values, even noise arrays can see some value from compression.

Encoding flat labels is fast.

### Condensed (Variable Width) Pins Label Format

`| BACKGROUND COLOR (STORED_DATA_WIDTH) | NUM_LABELS (u64) | UNIQUE LABELS (NUM_LABELS \* STORED_DATA_WIDTH) | FMT_BYTE | PIN SECTION |`

FMT_BYTE: `0000DDNN`

DD: 2^(DD) is the depth width
NN: 2^(NN) is the num pins width

PIN SECTION: `| PINS FOR LABEL 0 | PINS FOR LABEL 1 | ... | PINS FOR LABEL N |`

PINS: `| num_pins | INDEX_0 | INDEX_1 | ... | INDEX_N | DEPTH_0 | DEPTH_1 | ... | DEPTH_N |`

Note that INDEX_0 to INDEX_N are stored with a difference filter applied to improve compressibility.

A pin (color, position, depth) is a line segment that joins together multiple connected component IDs and labels them with a color (an index into UNIQUE LABELS) in order to use 3D information to compress the labels as compared with the flat label format. Pins are slow to compute but fast to decode, however random access is lost (a full scan of the labels section is needed to decode a subset of crack codes). The most frequent pin is replaced with a background color. Like with flat, efficient reading of unique labels, efficient remapping, and search are supported. 

Depending on the image statistics and quality of the pin solver, pins can be much smaller than flat or larger (some heuristics are used to avoid this case). An excellent example of where pins do well is a binary image where remarkable savings can be achieved in the labels section (though overall it is probably a small part of the file).

### Crack Code Format


CRACK CODE: `| CHAIN 0 | CHAIN 1 | ... | CHAIN N |`

CHAIN: `| BEGINNING OF CHAIN INDEX (sizeof(sx * sy)) | BIT PACKED MOVES (2 bits each) |`

The BEGINNING OF CHAIN INDEX (BOC) locates the grid vertex where the crack code will begin. Vertices are the corners of the pixel grid, with 0 at the top left and sx\*sy-1 at the bottom right (fortran order). 

The crack code is a NEWS code (up,right,left,down). Impossible combinations of directions are used to signal branching and branch termination. The next chain begins in the next byte when a termination signal causes the current branch count to reach zero.

There may be ways to further improve the design of the crack code. For example, by applying a difference filter a few more percent compression under gzip can be obtained. In the literature, there are other shorter codes such as a left,right,straight (LRS) code and fancy large context compressors that can achieve fewer than one bit per a move.

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

3. Freeman, H., 1974. Computer Processing of Line-Drawing Images. ACM Comput. Surv. 6, 57–97. https://doi.org/10.1145/356625.356627


