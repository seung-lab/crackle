from collections import defaultdict
from typing import List

import numpy as np
import networkx as nx
from tqdm import tqdm

from .lib import compute_dtype
from .ccl import color_connectivity_graph

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

def create_graph(
  labels, permissible=False, include_borders=False
):
  sx, sy = labels.shape
  G = nx.Graph()

  # graph is of corners and edges
  # origin is located at top left
  # corner of the image

  sxe = sx + 1 # sx edges
  sye = sy + 1 # sy edges

  check = lambda a,b: a != b
  if permissible:
    check = lambda a,b: a == b

  # assign vertical edges
  for x in range(1, sx):
    for y in range(0, sy):
      if check(labels[x,y], labels[x-1,y]):
        node_up = x + sxe * y
        node_down = x + sxe * (y + 1)
        G.add_edge(node_up, node_down)

  # assign horizontal edges
  for x in range(0, sx):
    for y in range(1, sy):
      if check(labels[x,y], labels[x,y-1]):
        node_left = x + sxe * y
        node_right = (x+1) + sxe * y
        G.add_edge(node_left, node_right)

  if include_borders and not permissible:
    for x in range(1,sxe): # vertical
      G.add_edge(x-1, x)
      G.add_edge(x-1 + sxe * (sye-1), x + sxe * (sye-1))

    for y in range(1,sye): # horizontal
      G.add_edge(sxe * y, sxe * (y+1))
      G.add_edge((sxe-1) + sxe * y, (sxe-1) + sxe * (y+1))

  return G

def create_crack_codes(labels, permissible) -> List[List[int]]:
  sx, sy = labels.shape
  G = create_graph(labels, permissible=permissible)
  Gcc = list(nx.connected_components(G))
  
  if len(Gcc) == 0:
    return []

  dirmap = {
    1: 'r',
    -1: 'l',
    (sx+1): 'd',
    -(sx+1): 'u',
  }

  revisit = []
  chains = []

  for i in range(len(Gcc)):
    cc = Gcc[i]
    remaining = set([ 
      tuple(sorted(edg)) for edg in G.subgraph(cc).edges 
    ])
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
        code.append('t')
        branches_taken -= 1
        if len(revisit):
          node = revisit.pop()
        elif len(remaining):
          node = next(iter(remaining))[0]
        continue
      elif len(neighbors) > 1:
        code.append('b')
        revisit.append(node)
        branches_taken += 1

      next_node = neighbors.pop()
      dir_taken = dirmap[next_node - node]
      code.append(dir_taken)
      remaining.discard(tuple(sorted([node,next_node])))
      node = next_node

    while branches_taken > 0:
      code.append('t')
      branches_taken -= 1

    chains.append(code)

  return symbols_to_integers(chains)

def symbols_to_integers(chains):
  encoded_chains = []

  LEFT = 0b10
  RIGHT = 0b01
  UP = 0b00
  DOWN = 0b11
  
  BRANCH = [UP,DOWN]
  BRANCH2 = [LEFT,RIGHT]
  TERM = [DOWN,UP]
  TERM2 = [RIGHT,LEFT]

  for chain in chains:
    code = [ chain[0] ] # node
    for symbol in chain[1:]:
      if symbol == 'u':
        code.append(UP)
      elif symbol == 'd':
        code.append(DOWN)
      elif symbol == 'l':
        code.append(LEFT)
      elif symbol == 'r':
        code.append(RIGHT)
      elif symbol == 'b':
        if code[-1] != BRANCH[1]:
          code.extend(BRANCH)
        else:
          code.extend(BRANCH2)
      elif symbol == 't':
        if code[-1] != TERM[1]:
          code.extend(TERM)
        else:
          code.extend(TERM2)
      else:
        raise ValueError(f"Got unsupported symbol: {symbol}")
    encoded_chains.append(code)

  return encoded_chains

def decode_crack_code(chains, sx, sy):
  # voxel connectivity
  # four bits: -y-x+y+x true is passable
  edges = np.full((sx,sy), fill_value=0b1111, dtype=np.uint8) 

  left = 0b10
  right = 0b01
  up = 0b00
  down = 0b11

  sxe = sx + 1
  sye = sy + 1

  # graph is of corners and edges
  # origin is located at top left
  # corner of the image
  for node, symbols in chains.items():
    y = node // sxe
    x = node - (sxe * y)

    revisit = []
    for symbol in symbols:
      if symbol == up:
        if (x-1) < sx and (y-1) < sy:
          edges[x-1,y-1] = edges[x-1,y-1] & 0b1110
        if x < sx and (y-1) < sy:
          edges[x,y-1] = edges[x,y-1] & 0b1101
        y -= 1
      elif symbol == down:
        if (x-1) < sx and y < sy:
          edges[x-1,y] = edges[x-1,y] & 0b1110
        if x < sx and y < sy:
          edges[x,y] = edges[x,y] & 0b1101
        y += 1
      elif symbol == left:
        if (x-1) < sx and (y-1) < sy:
          edges[x-1,y-1] = edges[x-1,y-1] & 0b1011
        if (x-1) < sx and y < sy:
          edges[x-1,y] = edges[x-1,y] & 0b0111
        x -= 1
      elif symbol == right:
        if x < sx and (y-1) < sy:
          edges[x,y-1] = edges[x,y-1] & 0b1011
        if x < sx and y < sy:
          edges[x,y] = edges[x,y] & 0b0111
        x += 1
      elif symbol == 'b':
        revisit.append((x,y))
      elif symbol =='t':
        if len(revisit) > 0:
          x, y = revisit.pop()

  return color_connectivity_graph(edges)

def unpack_binary(code, sx, sy):
  chains = defaultdict(list)

  if len(code) == 0:
    return chains

  dtype = compute_dtype((sx+1) * (sy+1))
  num_moves = np.dtype(dtype).itemsize * 8 // 2
  
  code = np.frombuffer(code, dtype=dtype)

  branches_taken = 0
  node = 0
  symbols = []
  for moveset in code:
    if branches_taken == 0:
      node = int(moveset)
      branches_taken = 1
      continue

    for i in range(num_moves):
      move = (moveset >> (2*i)) & 0b11
      if (
        (move == 0 and len(symbols) > 0 and symbols[-1] == 3)
        or (move == 2 and len(symbols) > 0 and symbols[-1] == 1)
      ):
        symbols[-1] = 't' # terminate
        branches_taken -= 1
        if branches_taken == 0:
          break
      elif (
        (move == 3 and len(symbols) > 0 and symbols[-1] == 0)
        or (move == 1 and len(symbols) > 0 and symbols[-1] == 2)
      ):
        symbols[-1] = 'b' # branch
        branches_taken += 1
      else:
        symbols.append(move)

    if branches_taken == 0:
      chains[node].extend(symbols)
      symbols = []

  return chains

def pack_codes(
  chains:List[List[int]], sx:int, sy:int
) -> bytes:
  binary = b''

  dtype = compute_dtype((sx+1) * (sy+1))
  num_bytes = np.dtype(dtype).itemsize
  num_moves = num_bytes * 8 // 2

  for chain in chains:
    node = dtype(chain.pop(0))
    binary += node.tobytes()[:num_bytes]
 
    for moveset in sip(chain, num_moves):
      encoded = dtype(0)
      for i, move in enumerate(moveset):
        encoded |= (move << (2*i))
      binary += encoded.tobytes()[:num_bytes]

  return binary

def encode_boundaries(labels, permissible:bool  =False):
  sx, sy, sz = labels.shape

  binary_components = []
  for z in range(sz):
    codes = create_crack_codes(labels[:,:,z], permissible)
    binary_components.append(
      pack_codes(codes, sx, sy)
    )

  return binary_components

