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

@click.command()
@click.option("-c/-d", "--compress/--decompress", default=True, is_flag=True, help="Compress from or decompress to a numpy .npy file.", show_default=True)
@click.option('-i', "--info", default=False, is_flag=True, help="Print the header for the file.", show_default=True)
@click.option('-l', "--labels", default=False, is_flag=True, help="Print unique labels contained in the image.", show_default=True)
@click.option('-t', "--test", default=False, is_flag=True, help="Check for file corruption and report damaged areas.", show_default=True)
@click.option('-p', '--allow-pins', default=False, is_flag=True, help="Allow pin encoding.", show_default=True)
@click.option('-m', '--markov', default=0, help="If >0, use this order of markov compression for the crack code.", show_default=True)
@click.option('-k', '--keep', default=False, is_flag=True, help="Keep the original file.", show_default=True)
@click.option('-z', 'gzip', default=False, is_flag=True, help="Apply gzip compression after encoding.", show_default=True)
@click.argument("source", nargs=-1)
def main(compress, info, test, labels, allow_pins, markov, source, keep, gzip):
	"""
	Compress and decompress crackle (.ckl) files to and from numpy (.npy) files.

	Compatible with crackle format version 0 streams.
	"""
	for i in range(len(source)):
		if source[i] == "-":
			source = source[:i] + sys.stdin.readlines() + source[i+1:]
	
	for src in source:
		if info:
			print_header(src)
			continue
		elif test:
			check_binary(src)
			continue
		elif labels:
			print_labels(src)
			continue

		if compress:
			compress_file(src, allow_pins, markov, gzip, keep)
		else:
			decompress_file(src, keep)

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

	print(f"Filename: {src}")
	for key,val in head.__dict__.items():
		print(f"{key}: {val}")
	print(f"num_labels: {num_labels}")
	print()

def decompress_file(src, keep):
	try:
		arr = crackle.util.aload(src)
	except FileNotFoundError:
		print(f"crackle: File \"{src}\" does not exist.")
		return
	except crackle.FormatError as err:
		print("crackle:", err)
		return
	except crackle.DecodeError:
		print(f"crackle: {src} could not be decoded.")
		return

	dest = src.replace(".ckl", "").replace(".gz", "").replace(".xz", "").replace(".lzma", "")
	_, ext = os.path.splitext(dest)
	
	if ext != ".npy":
		dest += ".npy"

	crackle.save_numpy(arr, dest)

	try:
		stat = os.stat(dest)
		if stat.st_size == 0:
			raise ValueError("File is zero length.")
		if not keep:
			os.remove(src)
	except (FileNotFoundError, ValueError) as err:
		print(f"crackle: Unable to write {dest}. Aborting.")
		sys.exit()

def compress_file(src, allow_pins, markov, gzip, keep):
	is_crackle = False
	try:
		data = crackle.util.load_numpy(src)
	except ValueError:
		try:
			data = crackle.aload(src, allow_mmap=True)
			is_crackle = True
		except:
			print(f"crackle: {src} is not a numpy or crackle file.")
			return
	except FileNotFoundError:
		print(f"crackle: File \"{src}\" does not exist.")
		return

	orig_src = src
	src = removesuffix(src, ".lzma")
	src = removesuffix(src, ".gz")
	src = removesuffix(src, ".xz")
	src = removesuffix(src, ".ckl")

	dest = f"{src}.ckl"
	if gzip:
		dest += ".gz"

	if is_crackle:
		data.binary = crackle.codec.reencode(data.binary, markov_model_order=int(markov))
		data.save(dest)
	else:
		crackle.save(data, dest, allow_pins=allow_pins, markov_model_order=int(markov))
	del data

	if dest == orig_src:
		return
		
	try:
		stat = os.stat(dest)
		if stat.st_size == 0:
			raise ValueError("File is zero length.")
		if not keep:
			os.remove(orig_src)
	except (FileNotFoundError, ValueError) as err:
		print(f"crackle: Unable to write {dest}. Aborting.")
		sys.exit()

def removesuffix(x:str, suffix:str) -> str:
  if x.endswith(suffix):
    x = x[:-len(suffix)]
  return x
