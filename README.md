[![PyPI version](https://badge.fury.io/py/crackle-codec.svg)](https://badge.fury.io/py/crackle-codec)

# Crackle: Next gen. 3D segmentation compression codec.

```bash
# Command Line Interface
crackle data.npy # creates data.ckl
crackle -m 5 data.npy # use a 5th order context model
crackle -p data.npy # use pin encoding for labels
crackle -p -m 5 data.npy # use pins and markov model
crackle -d data.ckl # recovers data.npy
crackle -m 0 data.ckl # change markov model order
crackle -t data.ckl # check for file corruption
```

```python
import crackle
import numpy

labels = np.load("example.npy") # a 2D or 3D dense segmentation

binary = crackle.compress(labels, allow_pins=False, markov_model_order=0)
labels = crackle.decompress(binary, parallel=0) # use all cores (default)

# faster extraction of binary images
binary_image = crackle.decompress(binary, label=1241)

# get unique labels without decompressing
uniq = crackle.labels(binary) 
# get num labels without decompressing
N = crackle.num_labels(binary) 
# get min and max without decompressing
mn = crackle.min(binary)
mx = crackle.max(binary)
# check if label in array in log(N) time
has_label = crackle.contains(binary, label)

# Remap labels without decompressing. Could
# be useful for e.g. proofreading.
remapped = crackle.remap(
  binary, { 1: 2, 2: 3, ... },
  preserve_missing_labels=True
)

# change dtype to smallest possible w/o precision loss
remapped = crackle.refit(binary)
# renumber array and change dtype to smallest possible
remapped = crackle.renumber(binary, start=0)

# for working with files
# if .gz is appended to the filename, the file will be
# automatically gzipped (or ungzipped)
crackle.save(labels, "example.ckl.gz")
labels = crackle.load("example.ckl.gz")

# Save a crackle array as a numpy array
# in a memory efficient manner.
crackle.save(binary, "example.npy.gz")

arr = crackle.CrackleArray(binary, parallel=0) # 0 means use all cores (default)
res = arr[:10,:10,:10] # array slicing (efficient z ranges)
arr[:,:,30] = 20 # write to a crackle array (whole z slices write faster)
20 in arr # log(N) check
arr = arr.numpy() # convert to a numpy array

# low memory extraction of point clouds
ptc = crackle.point_cloud(binary) # { label: np.ndarray, ... }
ptc = crackle.point_cloud(binary, label=777)
ptc = crackle.point_cloud(binary) 

# rapid and low memory
voxel_counts = crackle.voxel_counts(binary)

# building big arrays with low memory
binary = crackle.zeros([5000,5000,5000], dtype=np.uint64, order='F')

part1 = np.zeros([1000, 1000, 1000], dtype=np.uint32)
part2 = crackle.ones([1000, 1000, 1000], dtype=np.uint32)

binary = crackle.asfortranarray(binary)
binary = crackle.ascontiguousarray(binary)

# creates a crackle binary with part1 stacked atop part2
# in the z dimension. x and y dimensions must match
# without needing to decompress anything.
binary = crackle.zstack([ part1, part2 ])

# splits a crackle binary into before, middle (single slice),
# and after sections without decompressing.
before, middle, after = crackle.zsplit(binary, z=742)

# splits binary into individual z slices
sections = crackle.zshatter(binary)
```

*This repository is currently Beta. It works and the format is reasonably fixed. There may be some improvements down the line (such as 3d compression of crack codes), but they will be a new format version number.*

Crackle is a compression codec for 3D dense segmentation (labeled) images. The algorithm accepts both signed and unsigned integer labels (though the implementation currently has some restrictions on signed integers). It is written in C++ and has Python bindings. Crackle uses a two pass compression strategy where the output of crackle may be further comrpessed with a bitstream compressor like gzip, bzip2, zstd, or lzma. However, if the Crackle binary, which is already small, is not further compressed, it supports several efficient operations:

- Query if a label exists in the image
- Extract unique labels
- Remap labels
- Decode by Z-Range

Crackle is inspired by Compresso \[1\]. Compresso innovated by separating labels from boundary structures. There were conceptually four (but really five) elements in the format: header, labels, bit packed and RLE encoded binary image boundaries, and indeterminate boundary locations. 

Crackle improves upon Compresso by replacing the bit-packed boundary map with a "crack code" and can also use 3D information to reduce redundancy in labels using "pins".

See benchmarks for more information on Crackle's size and compute effiency.

## Installation

```bash
pip install crackle-codec 
```

Building from source (requires cmake and a c++ compiler):

```
git clone https://github.com/seung-lab/crackle.git
cd crackle
git submodule update --init --recursive # fetches google/crc32c library
python setup.py develop
```

## Versions

| Format Version | Description                                                    |
|----------------|----------------------------------------------------------------|
| 0              | Initial release w/ flat, pins, crack codes with finite context modeling. Beta. |
| 1              | Incr. header to 29 bytes from 24. num_label_bytes u32->u64, adds crcs to protect stream components. |

## Stream Format

| Section   | Bytes                                     | Description                                                                                                     |
|-----------|-------------------------------------------|-----------------------------------------------------------------------------------------------------------------|
| Header    | v0: 24, v1: 29                            | Metadata incl. length of fields.                                                                                |
| Crack Index     | header.sz * sizeof(uint32), v1: +4            | Number of bytes for the crack codes in each slice + CRC32c(le)
| Labels       | header.num_label_bytes        | Can be either "flat" labels or "pins". Describes how to color connected components.                                                                                   |
| Crack Codes    | Variable length.           | Instructions for drawing crack boundaries.             |
| Labels crc32c    | (v1 only) 4(le)           | v0: n/a, v1: crc32c of the labels binary.             |
| Labels crc32c    | (v1 only) header.sz * 4(le)           | v0: n/a, v1: crc32c of the uncompressed uint32_t fortran order connected component labeling of each z-slice.             |

### A Note on CRCs

CRCs protect each step of the decoding process. The fixed width header is protected by crc8, which contains information for decoding the crack index. The crack index is in turn protected by a crc32c. This is not overkill because a very large volume or a volume that is randomly accessible in XY as well as Z would need a crc32 vs a crc16.

The crack index is used for decoding the structural information (the connected components for each slice). We store a crc32c for each z-slice. This allows random z access to be validated while balancing against the storage cost of creating too many crc (e.g. vs. once per a grid).

We also store a crc32c for the labels binary independently of the crack codes.

All crcs are stored little endian.

Why not store a single crc for the entire uncompressed image? This would make it difficult to validate as a single crackle file could represent many terabytes of data. It would also make it difficult to edit the labels (remap) independently of the structure. Storing a crc32c per a z-slice also allows for z-stacking independent slices without recalculating crcs.

The downside to this strategy is a small increase in the file size and an increase in false positives for crc32s. This is the price of having the format be more of a random-access array format than a bitstream format. However, as crackle is designed to be two stage compressed, for example, with `lzip`, an LZMA variant with error correction properties, these issues are mitigated when archived.

crc8 (0xe7, initialized with 0xFF) was selected due to its ability to reliably detect [two bit flips](https://users.ece.cmu.edu/~koopman/crc/index.html) in up to 247 bits of message data, the best available for our header length.

crc32c was selected as the polynomial due to the availability of [high performance implementations](https://dougallj.wordpress.com/2022/05/22/faster-crc32-on-the-apple-m1/). This is important to avoid CRC calculation being a significant cost to the codec.

#### Error Detection and (Limited, Human Assisted) Correction

Due to this mutli-crc strategy, it is possible to narrow down corruptions to the section of the binary where they occur. For example, if you are concerned with only z=1-100 and the error occurs at z=200, you're ok. If the error occurs in the labels_binary, but you were planning on applying a full new mapping anyway, you can get away with discarding the extant labeling. Certain bit flips in the labels binary will create out of range keys, which will aid in identifying exactly where the error occured.  Headers can be repaired potentially be human inspection (if they know the dataset).

### Header


| Attribute         | Value             | Type    | Description                                     |
|-------------------|-------------------|---------|-------------------------------------------------|
| magic             | crkl              | char[4] | File magic number.                              |
| format_version    | 0 or 1            | u8      | Stream version.                   |
| format_field      | bitfield          | u16     | See below.                 |
| sx, sy, sz        | >= 0              | u32 x 3 | Size of array dimensions.                       |
| grid_size         | log2(grid_size)   | u8      | Stores log2 of grid dimensions in voxels.          |
| num_label_bytes   | Any.              | u64      | Number of bytes of the labels section. Note the labels come in at least two format types.          |
| crc8             | Any.              | u8      | CRC8 of format_field thru num_label_bytes using polynomial 0xe7 (implicit) and 0xFF initialization. |


Format Field (u16): DDSSCLLFGOOOOURR (each letter represents a bit, left is LSB)

DD: 2^(DD) = byte width of returned array (1,2,4,8 bytes)  
SS: 2^(SS) = byte width of stored labels (sometimes you can store values in 2 bytes when the final array is 8 bytes)  
C: 1: crack codes denote impermissible boundaries 0: they denote permissible boundaries.  
LL: 0: "flat" label format, 1: fixed width pins (unused?) 2: variable width pins 3: reserved  
F: whether the array is to be rendered as C (0) or F (1) order  
G: Signed (if (1), data are signed int, otherwise unsigned int)  
OOOO: Nth-Order of Markov Chain (as an unsigned integer, typical values 0, or 3 to 7). If 0, markov compression is disabled.  
U: if 0, unique labels are sorted, else, unsorted  
R: Reserved  

CRC8 only covers the header. It doesn't cover the magic number or format version since those are easily human correctable if needed.

### Flat Label Format

| Attribute     | Type                                        | Description                                                                                                 |
|---------------|---------------------------------------------|-------------------------------------------------------------------------------------------------------------|
| num_unique    | u64                                         | Number of unique labels in this volume.                                                                     |
| unique_labels | stored_type[num_unique]                     | Sorted ascending array of all unique values in image, stored in the smallest data type that will hold them. |
| cc_per_grid   | smallest_type(sx \* sy)[sz]                 | Array containing the number of CCL IDs in each grid (usually a z-slice).                                    |
| cc_to_labels  | smallest_type(num_labels)[sum(cc_per_grid)] | Array mapping CCL IDs to their proper value by indexing the unique labels array.                            |

Flat labels are random access read, allow efficient reading of unique labels, efficient remapping, and efficient search for a given label's existence. Since the connected component labels can often use a smaller byte width than the unique values, even noise arrays can see some value from compression.

Encoding flat labels is fast.

### Condensed (Variable Width) Pins Label Format

| Attribute        | Type                                | Description                                                                                                 |
|------------------|-------------------------------------|-------------------------------------------------------------------------------------------------------------|
| background_color | stored_data_width                   | Background color of image.                                                                                  |
| num_unique       | u64                                 | Number of unique labels in this volume.                                                                     |
| unique_labels    | stored_type[num_unique]             | Sorted ascending array of all unique values in image, stored in the smallest data type that will hold them. |
| cc_per_grid   | smallest_type(sx \* sy)[sz]                 | Array containing the number of CCL IDs in each grid (usually a z-slice).                                    |
| fmt_byte         | u8                                  | 00CCDDNN  DD: 2^(DD) is the depth width NN: 2^(NN) is the num pins width, CC: 2^(CC) is the single components width.  |
| pin_section      | Bitstream to end of labels section. | Contains pin information.                                                                                   |

PIN SECTION: `| PINS FOR LABEL 0 | PINS FOR LABEL 1 | ... | PINS FOR LABEL N |`

PINS: `| num_pins | INDEX_0 | INDEX_1 | ... | INDEX_N | DEPTH_0 | DEPTH_1 | ... | DEPTH_N | num_single_labels | CC 0 | CC 1 | ... | CC N |`

Both `num_pins` and `num_single_labels` use the `num_pins_width`.

Note that INDEX_0 to INDEX_N are stored with a difference filter applied to improve compressibility.

A pin (color, position, depth) is a line segment that joins together multiple connected component IDs and labels them with a color (an index into UNIQUE LABELS) in order to use 3D information to compress the labels as compared with the flat label format. Pins are slow to compute but fast to decode, however random access is lost (a full scan of the labels section is needed to decode a subset of crack codes). The most frequent pin is replaced with a background color. Like with flat, efficient reading of unique labels, efficient remapping, and search are supported. 

Depending on the image statistics and quality of the pin solver, pins can be much smaller than flat or larger (some heuristics are used to avoid this case). An excellent example of where pins do well is a binary image where remarkable savings can be achieved in the labels section (though overall it is probably a small part of the file).

For very short pins (e.g. depth 0 or 1) that take more bytes to record than simply listing the corresponding CC label, we list the CC label instead. This calculation is made depending on the dimensions of the image and the max pin depth, and the byte width of the CCL labels.

Example calculation. For a 512 x 512 x 32 file with an average of 1000 CCL's per a slice and a maximum pin depth of 30, a pin takes 4 index + 1 depth = 5 bytes while a CCL takes 2 bytes. Therefore, depth 1 and 2 pins can be efficiently replaced with 1 and 2 CCL labels for a 60\% and 20\% savings respectively. CCLs are also difference coded to enhance second stage compressibility.

### Fixed Width Pins (disabled)

`| BACKGROUND COLOR (STORED_DATA_WIDTH) | NUM_LABELS (u64) | UNIQUE LABELS (NUM_LABELS \* STORED_DATA_WIDTH) | PIN SECTION |`

PIN SECTION: `|PIN0|PIN1|PIN2|...|PINN|`
PIN: `|LABEL|INDEX|DEPTH|`

A fixed width variant of pins has also been developed but is not enabled. It frequently is not significantly smaller than flat outside of special circumstances such as a binary image. An advantage this format would have over condensed is that the pins can be sorted and searched rapidly by index, which reduces the amount of reading one might have to do on an mmapped file. Please raise an issue if this seems like something that might be useful to you.

### Crack Code Format

CRACK CODE: `MARKOV MODEL | CHAIN 0 | CHAIN 1 | ... | CHAIN N |`

CHAIN: `| BEGINNING OF CHAIN INDEX (sizeof(sx * sy)) | BIT PACKED MOVES (2 bits each) |`

MARKOV MODEL (if enabled): `priority order of moves UDLR packed per a byte`. 4^order bytes.

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


