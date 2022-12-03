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

  for x in range(1, sx):
    for y in range(0, sy):
      if labels[x,y] != labels[x-1,y]:
        node_up = x + sx * y
        node_down = x + sx * (y + 1)
        G.add_edge(node_up, node_down)

  for x in range(0, sx):
    for y in range(1, sy):
      if labels[x,y] != labels[x,y-1]:
        node_left = x + sx * y
        node_right = (x+1) + sx * y
        G.add_edge(node_left, node_right)

  if include_borders:
    for x in range(1,sx):
      G.add_edge(x-1, x)
      G.add_edge(x-1 + sx * (sy-1), x + sx * (sy-1))

    for y in range(1,sy):
      G.add_edge(sx * y, sx * (y+1))
      G.add_edge((sx-1) + sx * y, (sx-1) + sx * (y+1))

  return G

def create_crack_code(labels):
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
    sx: down,
    -sx: up,
  }

  # special codes
  BRANCH = [up,down]
  TERM = [down,up]
  UNUSED1 = [left,right]
  UNUSED2 = [right,left]

  visited = set()
  revisit = []

  chains = []

  node = list(Gcc[0])[0]
  for i in range(len(Gcc)):
    cc = Gcc[i]
    remaining = set(cc)
    node = next(iter(remaining))
    remaining.discard(node)

    code = [ node ]
    while len(remaining) or len(revisit):
      neighbors = [ 
        n for n in G.neighbors(node) if n not in visited 
      ]
      visited.add(node)
      remaining.discard(node)

      if len(neighbors) == 0:
        code.extend(TERM)
        if len(revisit):
          node = revisit.pop()
        elif len(remaining):
          node = next(iter(remaining))
        continue
      elif len(neighbors) > 1:
        code.extend(BRANCH)
        revisit.append(node)

      next_node = neighbors.pop()
      dir_taken = dirmap[next_node - node]
      code.append(dir_taken)
      node = next_node
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
    x = node - sx * y

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

def unpack_crack_binary(code):
  symbols = []

  if len(chain) == 0:
    return symbols

  code = np.frombuffer(code, dtype=np.uint32)
  chains = {}

  branches_taken = 0
  node = 0
  for moveset in chain:
    if branches_taken == 0:
      node = int(moveset)
      branches_taken = 1
      continue

    symbols = []
    for i in range(16):
      move = (moveset >> (2*i)) & 0b11
      
      if move == 0 and symbols[-1] == 3 and len(symbols) > 1:
        symbols[-1] = 't' # terminate
        branches_taken -= 1
        if branches_taken == 0:
          break
      elif move == 3 and symbols[-1] == 0 and len(symbols) > 1:
        symbols[-1] = 'b' # branch
        branches_taken += 1
      else:
        symbols.append(move)
    chains[node] = symbols

  return chains

def crack_code_to_binary(chains):
  binary = b''

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
    chains = create_crack_code(labels[:,:,z])
    binary_components.append(
      crack_code_to_binary(chains)
    )

  return binary_components

