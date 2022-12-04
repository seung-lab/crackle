from collections import defaultdict
from typing import List

import numpy as np
import networkx as nx
from tqdm import tqdm

from .ccl import connected_components

def sip(iterable, block_size):
  """Sips a fixed size from the iterable."""
  ct = 0
  block = []
  for x in iterable:
    ct += 1
    block.append(x)
    if ct == block_size:
      yield block
      ct = 0
      block = []

  if len(block) > 0:
    yield block

def create_graph(labels, include_borders=False):
  sx, sy = labels.shape
  G = nx.Graph()

  # graph is of corners and edges
  # origin is located at top left
  # corner of the image

  sxe = sx + 1 # sx edges
  sye = sy + 1 # sy edges

  # assign vertical edges
  for x in range(1, sx):
    for y in range(0, sy):
      if labels[x,y] != labels[x-1,y]:
        node_up = x + sxe * y
        node_down = x + sxe * (y + 1)
        G.add_edge(node_up, node_down)

  # assign horizontal edges
  for x in range(0, sx):
    for y in range(1, sy):
      if labels[x,y] != labels[x,y-1]:
        node_left = x + sxe * y
        node_right = (x+1) + sxe * y
        G.add_edge(node_left, node_right)

  if include_borders:
    for x in range(1,sxe): # vertical
      G.add_edge(x-1, x)
      G.add_edge(x-1 + sxe * (sye-1), x + sxe * (sye-1))

    for y in range(1,sye): # horizontal
      G.add_edge(sxe * y, sxe * (y+1))
      G.add_edge((sxe-1) + sxe * y, (sxe-1) + sxe * (y+1))

  return G

def create_crack_codes(labels) -> List[List[int]]:
  sx, sy = labels.shape
  G = create_graph(labels)
  Gcc = list(nx.connected_components(G))
  
  if len(Gcc) == 0:
    return []

  left = 0b10
  right = 0b01
  up = 0b00
  down = 0b11

  dirmap = {
    1: right,
    -1: left,
    (sx+1): down,
    -(sx+1): up,
  }

  # special codes
  BRANCH = [up,down]
  TERM = [down,up]
  UNUSED1 = [left,right]
  UNUSED2 = [right,left]

  revisit = []
  chains = []

  for i in range(len(Gcc)):
    cc = Gcc[i]
    remaining = set(G.subgraph(cc).edges)
    node = next(iter(remaining))[0]
    remaining.discard(node)

    branches_taken = 1

    code = [ node ]
    while len(remaining) or len(revisit):
      neighbors = [ 
        n for n in G.neighbors(node) 
        if tuple(sorted([n,node])) in remaining
      ]
      
      if len(neighbors) == 0:
        code.extend(TERM)
        branches_taken -= 1
        if len(revisit):
          node = revisit.pop()
        elif len(remaining):
          node = next(iter(remaining))
        continue
      elif len(neighbors) > 1:
        code.extend(BRANCH)
        revisit.append(node)
        branches_taken += 1

      next_node = neighbors.pop()
      dir_taken = dirmap[next_node - node]
      code.append(dir_taken)
      remaining.discard(tuple(sorted([node,next_node])))
      node = next_node

    while branches_taken > 0:
      code.extend(TERM)
      branches_taken -= 1

    chains.append(code)

  return chains

def decode_crack_code(chains, sx, sy):
  # voxel connectivity
  # four bits: -y-x+y+x true is passable
  edges = np.zeros((sx,sy), dtype=np.uint8) 
  edges += 0b1111

  left = 0b10
  right = 0b01
  up = 0b00
  down = 0b11

  # graph is of corners and edges
  # origin is located at top left
  # corner of the image
  for node, symbols in chains.items():
    y = node // sx
    x = node - (sx * y)
    print(node,x,y)
    revisit = []
    for symbol in symbols:
      if symbol == up:
        edges[x-1,y-1] = edges[x-1,y-1] & 0b1110
        edges[x,y-1] = edges[x,y-1] & 0b1101
        y -= 1
      elif symbol == down:
        edges[x-1,y] = edges[x-1,y] & 0b1110
        edges[x,y] = edges[x,y] & 0b1101
        y += 1
      elif symbol == left:
        edges[x-1,y-1] = edges[x-1,y-1] & 0b1011
        edges[x-1,y] = edges[x-1,y] & 0b0111
        x -= 1
      elif symbol == right:
        edges[x,y-1] = edges[x,y-1] & 0b1011
        edges[x,y] = edges[x,y] & 0b0111
        x += 1
      elif symbol == 'b':
        revisit.append((x,y))
      elif symbol =='t':
        if len(revisit) > 0:
          x, y = revisit.pop()

  return connected_components(edges)

def unpack_binary(code):
  chains = defaultdict(list)

  if len(code) == 0:
    return chains

  code = np.frombuffer(code, dtype=np.uint32)

  branches_taken = 0
  node = 0
  for moveset in code:
    if branches_taken == 0:
      node = int(moveset)
      branches_taken = 1
      continue

    symbols = []
    for i in range(16):
      move = (moveset >> (2*i)) & 0b11
      if move == 0 and len(symbols) > 1 and symbols[-1] == 3:
        symbols[-1] = 't' # terminate
        branches_taken -= 1
        if branches_taken == 0:
          break
      elif move == 3 and len(symbols) > 1 and symbols[-1] == 0:
        symbols[-1] = 'b' # branch
        branches_taken += 1
      else:
        symbols.append(move)
    chains[node].extend(symbols)

  # print(chains)

  # for node, symbols in chains.items():
    # print(node, symbols)

  return chains

def pack_codes(chains:List[List[int]]) -> bytes:
  binary = b''

  # print("here", chains)

  for chain in chains:
    node = np.uint32(chain.pop(0))
    binary += node.tobytes()[:4]
 
    for moveset in sip(chain, 16):
      encoded = np.uint32(0)
      for i, move in enumerate(moveset):
        encoded |= (move << (2*i))
      binary += encoded.tobytes()[:4]

  return binary

def encode_boundaries(labels):
  sz = labels.shape[2]

  binary_components = []
  for z in tqdm(range(sz), desc='crack code z'):
    codes = create_crack_codes(labels[:,:,z])
    binary_components.append(
      pack_codes(codes)
    )

  return binary_components

