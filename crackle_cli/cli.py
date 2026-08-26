import os
import sys

import click
import numpy as np

import crackle

class Tuple3(click.ParamType):
  """A command line option type consisting of 3 comma-separated integers."""
  name = 'tuple3'
  def convert(self, value, param, ctx):
    if isinstance(value, str):
      try:
        value = tuple(map(int, value.split(',')))
      except ValueError:
        self.fail(f"'{value}' does not contain a comma delimited list of 3 integers.")
      if len(value) != 3:
        self.fail(f"'{value}' does not contain a comma delimited list of 3 integers.")
    return value

def normalize_extension(path):
	path = removesuffix(path, ".lzma")
	path = removesuffix(path, ".gz")
	path = removesuffix(path, ".xz")
	path = removesuffix(path, ".bz2")
	path = removesuffix(path, ".zstd")
	return path

def normalize_source(source):
	for i in range(len(source)):
		if source[i] == "-":
			source = source[:i] + sys.stdin.readlines() + source[i+1:]
	return source

@click.command()
@click.option('--allow-pins', default=False, is_flag=True, help="Allow pin encoding.", show_default=True)
@click.option('-m', '--markov', default=0, help="If >0, use this order of markov compression for the crack code.", show_default=True)
@click.option('-k', '--keep', default=False, is_flag=True, help="Keep the original file.", show_default=True)
@click.option('-z', 'gzip', default=False, is_flag=True, help="Apply gzip compression after encoding.", show_default=True)
@click.option('-p', '--parallel', default=0, help="Number of threads to use. 0 = num cpu", show_default=True)
@click.argument("source", nargs=-1)
def cckl(
	allow_pins, markov, source, 
	keep, gzip, parallel,
):
	"""
	Compress to crackle (.ckl) files 
	from other lossless segmentation file types.

	Supports: ckl, cpso, npy, nrrd, tiff

	Compatible with crackle format version 0 and 1 streams.

	See also dckl (decompress), crackle (all operations)
	"""
	source = normalize_source(source)

	for src in source:
		ckl_path = convert_file(src, allow_pins, markov, gzip, keep, "ckl", parallel)

@click.command()
@click.option('-k', '--keep', default=False, is_flag=True, help="Keep the original file.", show_default=True)
@click.option('-z', 'gzip', default=False, is_flag=True, help="Apply gzip compression after encoding.", show_default=True)
@click.option('--to', 'to_format', default="npy", type=str, help="Specify destination format. Supports: npy, nii, nrrd", show_default=True)
@click.option('-p', '--parallel', default=0, help="Number of threads to use. 0 = num cpu", show_default=True)
@click.argument("source", nargs=-1)
def dckl(
	source, 
	keep, gzip, parallel,
	to_format,
):
	"""
	Decompress crackle (.ckl) files 
	to other lossless segmentation file types.

	Supports: cpso, npy, nrrd, tiff

	Compatible with crackle format version 0 and 1 streams.

	See also cckl (compress), crackle (all operations)
	"""
	source = normalize_source(source)

	for src in source:
		ckl_path = convert_file(src, False, 0, gzip, keep, to_format, parallel)

@click.group()
def main():
	"""Perform operations on crackle files.

	Crackle files are a lossless compression format
	for multilabel 2D or 3D image segmentations. 

	See also: cckl (compress), dckl (decompress)
	"""
	pass

@main.command("info")
@click.argument("source", nargs=-1)
def infocmd(source):
	"""Print the header for the file."""

	source = normalize_source(source)

	for src in source:
		print_header(src)

@main.command("labels")
@click.argument("source", nargs=-1)
def labelscmd(source):
	"""Print unique labels contained in the image."""

	source = normalize_source(source)

	for src in source:
		print_labels(src)

@main.command("convert")
@click.option('--allow-pins', default=False, is_flag=True, help="Allow pin encoding.", show_default=True)
@click.option('-m', '--markov', default=None, help="If >0, use this order of markov compression for the crack code.", show_default=True)
@click.option('-k', '--keep', default=False, is_flag=True, help="Keep the original file.", show_default=True)
@click.option('-z', 'gzip', default=False, is_flag=True, help="Apply gzip compression after encoding.", show_default=True)
@click.option('-t', '--to', 'to_format', default="ckl", type=str, help="Specify destination format. Supports: npy, nii, nrrd, ckl", show_default=True)
@click.option('-p', '--parallel', default=0, help="Number of threads to use. 0 = num cpu", show_default=True)
@click.argument("source", nargs=-1)
def convertcmd(
	source, 
	markov, keep, gzip, to_format,
	parallel,
):
	"""
	Convert between crackle (.ckl) files and other file types.

	Supports: ckl, cpso, npy, nrrd, tiff

	Compatible with crackle format version 0 and 1 streams.
	"""
	source = normalize_source(source)

	for src in source:
		convert_file(src, allow_pins, markov, gzip, keep, to_format, parallel)


@main.command("test")
@click.argument("source", nargs=-1)
def testcmd(source):
	"""Check for file corruption and report damaged areas."""
	source = normalize_source(source)

	for src in source:
		check_binary(src)

@main.command("meta")
@click.argument("source", nargs=-1)
@click.option('-s', '--skip-empty', default=False, is_flag=True, help="Skip sidecar generation for empty (background only) crackle files.", show_default=True)
def metacmd(source, skip_empty):
	"""Create sidecar parquet files.

Creates .meta.parquet files with voxel counts, bounding boxes.
These files can be manually used for faster operations but are entirely optional.
"""
	source = normalize_source(source)

	for src in source:
		create_sidecar_file(src, skip_empty)

@main.command("stack")
@click.argument("source")
@click.option('-o', '--output', default="merged.ckl", help="File path to output to.", show_default=True)
def stackcmd(source, output):
	"""zstack all crackle files in the directory in alphabetical order."""
	stack_dir(source, output)

@main.command("downsample")
@click.argument("source")
@click.option('-o', '--output', default="downsampled.ckl", help="File path to output to.", show_default=True)
@click.option('-s', '--sparse', is_flag=True, default=False, help="Don't count background in the mode.", show_default=True)
@click.option('-p', '--parallel', default=0, help="Number of threads to use. 0 = num cpu", show_default=True)
@click.option('-n', 'N', default=1, help="Number of times to serially downsample.", show_default=True)
def downsamplecmd(source, output, sparse, parallel, N):
	"""Use 2x2x1 mode pooling to downsample this image."""
	arr = crackle.aload(source, allow_mmap=True)
	arr.parallel = parallel

	arr_ds = arr
	for i in range(N):
		arr_ds = arr.mode_pooling_2x2x1(sparse=sparse)
		arr_ds.parallel = parallel

	arr_ds.save(output)

@main.command("view")
@click.argument("source")
@click.option('-z', default=None, type=str, help="Which z-slices to visualize.", show_default=True)
def viewcmd(source, z):
	"""View this image."""
	slc = (slice(None), slice(None), slice(None))
	if z is not None:
		if ',' in z:
			zstart, zend = z.split(',')
			zstart = int(zstart)
			zend = int(zend)
		else:
			zstart = int(z)
			zend = zstart + 1

		slc = (slice(None), slice(None), slice(zstart, zend))

	import microviewer
	arr = crackle.aload(source, allow_mmap=True)

	labels = arr[slc]

	microviewer.view(labels, seg=True)
	
def check_binary(src):
	try:
		arr = crackle.aload(src, allow_mmap=True)
	except FileNotFoundError:
		print(f"crackle: File \"{src}\" does not exist.")
		return
	except crackle.FormatError as err:
		print("crackle:", err)
		return

	print(f"testing {src}...")

	report = crackle.codec.check(arr.binary)

	def pretty(human, key):
		if report[key] == True:
			print(f"{human} ok.")
		elif report[key] == False:
			print(f"{human} damaged (or false positive crc check).")
		elif report[key] is None:
			print(f"{human} maybe ok. (no crc check in this format version)")

	pretty("header", "header")
	pretty("crack index", "crack_index")
	pretty("labelling", "labels")

	if report["z"] is None:
		print("sections maybe ok. (no crc check in this format version)")
	elif report["z"] == []:
		print("sections ok.")
	else:
		print(f"sections damaged: { ','.join(report['z']) }")

	print("done.")

def print_labels(src):
	try:
		arr = crackle.aload(src)
	except FileNotFoundError:
		print(f"crackle: File \"{src}\" does not exist.")
		return
	except crackle.FormatError as err:
		print("crackle:", err)
		return

	labels = arr.labels()
	print("\n".join(( str(l) for l in labels )))

def toSI(value) -> str:
	if value > 1e12:
		return f"{value / 1e12:.2f} tera"
	elif value > 1e9:
		return f"{value / 1e9:.2f} giga"
	elif value > 1e6:
		return f"{value / 1e6:.2f} mega"
	elif value > 1e3:
		return f"{value / 1e3:.2f} kilo"
	else:
		return f"{value} "

def print_header(src):
	try:
		head = crackle.util.load_header(src, ignore_crc_check=True)
	except FileNotFoundError:
		print(f"crackle: File \"{src}\" does not exist.")
		return
	except crackle.FormatError as err:
		print("crackle:", err)
		return

	num_labels = crackle.util.load_num_labels(src)
	magnitude = toSI(head.voxels()) + "voxels"

	num_total_bytes = toSI(os.path.getsize(src)) + "bytes"

	print(f"Filename: {src}")
	print(f"size: {magnitude}")
	print(f"bytes: {num_total_bytes}")
	for key,val in head.__dict__.items():
		print(f"{key}: {val}")
	print(f"num_labels: {num_labels}")
	print()

def create_sidecar_file(src, skip_empty):
	arr = crackle.aload(src, allow_mmap=True)

	if skip_empty and (arr.num_labels() <= 1) and (0 in arr):
		return

	meta_path = normalize_extension(src)
	meta_path += ".meta.parquet"
	arr.cache_meta(meta_path)

def convert_file(
	src,
	allow_pins:bool,
	markov:int,
	gzip:bool,
	keep:bool,
	to_format:str,
	parallel:int = 0,
):
	nsrc = normalize_extension(src)
	ckl_path = None

	try:
		if nsrc.endswith(".ckl"):
			arr = crackle.util.aload(src)
		else:
			arr = crackle.util.load_any(src)
	except FileNotFoundError:
		print(f"crackle: File \"{src}\" does not exist.")
		return None
	except crackle.FormatError as err:
		print("crackle:", err)
		return None
	except crackle.DecodeError:
		print(f"crackle: {src} could not be decoded.")
		return None
	except Exception as err:
		print(f"crackle: {src} could not be decoded.")
		print(err)
		return None

	dest, ext = os.path.splitext(nsrc)

	if to_format == "ckl":
		dest += ".ckl"
		if gzip:
			dest += ".gz"
		if isinstance(arr, crackle.CrackleArray):
			if arr.header().format_version == 0:
				# This causes file corruption for some reason. Requires more extensive inquiry.
				print("Reencoding format version 0 is not currently supported.")
				return

			arr.binary = crackle.codec.reencode(
				arr.binary,
				markov_model_order=int(markov),
				parallel=parallel,
			)
			arr.save(dest)
		else:
			crackle.save(
				arr, dest, 
				allow_pins=allow_pins, 
				markov_model_order=int(markov),
				parallel=parallel,
			)
		ckl_path = dest
	elif to_format == "npy":
		if not dest.endswith(".npy"):
			dest += ".npy"
		if gzip:
			dest += ".gz"
		crackle.save_numpy(arr, dest)
	elif to_format == "nrrd":
		if not dest.endswith(".nrrd"):
			dest += ".nrrd"
		compress = "gzip" if gzip else "raw"
		crackle.util.save_nrrd(arr, dest, compress=compress)
	elif to_format in ("tiff", "tif"):
		dest += f".{to_format}"
		compress = "zlib" if gzip else None
		crackle.util.save_tiff(arr, dest, compression=compress)
	elif to_format == "cpso":
		if not dest.endswith(".cpso"):
			dest += f".cpso"
		if gzip:
			dest += ".gz"
		crackle.util.save_compresso(arr, dest)

	if src == dest:
		keep = True

	if ckl_path is None and keep and to_format == "ckl":
		ckl_path = src

	try:
		stat = os.stat(dest)
		if stat.st_size == 0:
			raise ValueError("File is zero length.")
		if not keep:
			os.remove(src)
	except (FileNotFoundError, ValueError) as err:
		print(f"crackle: Unable to write {dest}. Aborting.")
		sys.exit()

	return ckl_path

def removesuffix(x:str, suffix:str) -> str:
  if x.endswith(suffix):
    x = x[:-len(suffix)]
  return x

def stack_dir(src:str, output:str):
	filenames = [ 
		filename 
		for filename in os.listdir(src)
		if filename.endswith(".ckl") 
	]
	filenames.sort()

	binaries = [ 
		crackle.bload(filename)
		for filename in filenames 
	]

	fused_binary = crackle.zstack(binaries)
	del binaries

	with open(os.path.join(src, output), "wb") as f:
		f.write(fused_binary)



