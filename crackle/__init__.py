"""A 3D dense labeled image compressor.

Crackle uses a strategy similar to Freeman 
crack codes to encode the boundaries of 
large dense 3D segmentation images.

If the number of voxel pairs in the image
are fewer than half the size of the image, 
the cracks are encoded as permissible edges, 
if they are more than half, they are coded as 
impermissible edges. This makes crackle work
well on both highly structured and noisy 
images.

The crack code is a NSEW code requiring 2 bits
per a move with explicit branch and termination
symbols encoded using impossible direction pairs 
(left,right), (up,down). Beginning of chain is 
encoded as an integer index to the start position
on the crack grid.

The crack codes can be futher compressed using
a finite context model.

Labels are encoded as a set of unique labels
in sorted order followed by an index of connected
component ids into the unique labels. This gives
good performance on both structured and random 
images as often bytes needed to encode the connected
components are much fewer than the labels themselves.

There is also an alternative label encoding called
"pins" that draws a line segment connecting several
connected components at once that works well for
simpler images, but can be more expensive
on more complicated images. The pin encoding method is
experimental and its format is subject to change and
is by default disabled.

Author: William Silversmith
Affiliation: Princeton Neuroscience Institute
Date: December 2022 - March 2023
"""
from .array import CrackleArray, CrackleRemoteArray
from .codec import (
	compress, decompress, labels, remap, 
	nbytes, components, component_lengths,
	header, contains, crack_codes, refit,
	renumber, num_labels, min, max
)
from .headers import FormatError, CrackleHeader
from .util import save, load, aload

