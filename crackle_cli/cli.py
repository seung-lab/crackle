import os
import os.path
import sys

import click
import crackle
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
@click.option('--allow-pins', default=False, is_flag=True, help="Allow pin encoding.", show_default=True)
@click.option('-m', '--markov', default=0, help="If >0, use this order of markov compression for the crack code.", show_default=True)
@click.option('-z', 'gzip', default=False, is_flag=True, help="Apply gzip compression after encoding.", show_default=True)
@click.argument("source", nargs=-1)
def main(compress, info, allow_pins, markov, source, gzip):
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

		if compress:
			compress_file(src, allow_pins, markov, gzip)
		else:
			decompress_file(src)

def print_header(src):
	try:
		with open(src, "rb") as f:
			binary = f.read()
	except FileNotFoundError:
		print(f"crackle: File \"{src}\" does not exist.")
		return

	head = crackle.header(binary)
	print(f"Filename: {src}")
	for key,val in head.__dict__.items():
		print(f"{key}: {val}")
	print()

def decompress_file(src):
	try:
		data = crackle.util.load(src)
	except FileNotFoundError:
		print(f"crackle: File \"{src}\" does not exist.")
		return
	except crackle.DecodeError:
		print(f"crackle: {src} could not be decoded.")
		return

	dest = src.replace(".ckl", "").replace(".gz", "").replace(".xz", "").replace(".lzma", "")
	_, ext = os.path.splitext(dest)
	
	if ext != ".npy":
		dest += ".npy"

	np.save(dest, data)

	try:
		stat = os.stat(dest)
		if stat.st_size > 0:
			os.remove(src)
		else:
			raise ValueError("File is zero length.")
	except (FileNotFoundError, ValueError) as err:
		print(f"crackle: Unable to write {dest}. Aborting.")
		sys.exit()

def compress_file(src, allow_pins, markov, gzip):
	try:
		data = np.load(src)
	except ValueError:
		print(f"crackle: {src} is not a numpy file.")
		return
	except FileNotFoundError:
		print(f"crackle: File \"{src}\" does not exist.")
		return

	dest = f"{src}.ckl"
	if gzip:
		dest += ".gz"
	crackle.save(data, dest, allow_pins=allow_pins, markov_model_order=int(markov))
	del data

	try:
		stat = os.stat(dest)
		if stat.st_size > 0:
			os.remove(src)
		else:
			raise ValueError("File is zero length.")
	except (FileNotFoundError, ValueError) as err:
		print(f"crackle: Unable to write {dest}. Aborting.")
		sys.exit()

