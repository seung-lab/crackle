import numpy as np
import networkx as nx
from tqdm import tqdm

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

