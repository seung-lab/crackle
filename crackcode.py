import numpy as np
import compresso
import time
import networkx as nx
from cloudvolume import view
from cloudvolume.lib import sip
from tqdm import tqdm

import time

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
  sx,sy = labels.shape
  G = create_graph(labels)
  Gcc = list(nx.connected_components(G))
  
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

def estimate_slice(labels):
  sx, sy = labels.shape

  G = create_graph(labels)

  s = time.time()
  ncc = len(list(nx.connected_components(G)))
  # print(time.time() - s)

  n_edges = len(G.edges())
  ideal_size = 6 * ncc + (2 * n_edges) / 8
  cur_size = sx * sy / 8 / 2

  # print("best size: ", ideal_size, "bytes")
  # print("current est. size", cur_size, "bytes")
  # print("improvement", cur_size / ideal_size)
  return ideal_size

def estimate_improvement(labels):
  sx, sy, sz = all_labels.shape

  chain_size = 0
  for z in tqdm(range(sz)):
    chain_size += estimate_slice(all_labels[:,:,z])

  print(chain_size)

def pack(chains):
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

# def raw_labels(buf):
#   info = compresso.header(buf)

#   offset = 36#compresso.COMPRESSO_HEADER_SIZE
#   id_bytes = info["id_size"] * info["data_width"]
#   ldtype = compresso.label_dtype(info)
#   wdtype = compresso.window_dtype(info)

#   ids = np.frombuffer(buf[offset:offset+id_bytes], dtype=ldtype)
#   # print(ids.size)
#   return ids

with open("connectomics.npy.cpso", "rb") as f:
  cpso = f.read()
  all_labels = compresso.decompress(cpso)

total = 0
binary = b''
for z in tqdm(range(512)):
  chains = create_crack_code(all_labels[:,:,z])
  packed = pack(chains)
  binary += packed
  total += len(packed)
  # print(len(packed))
  # total += sum([ (4.0 if x > 100 else 0.25) for x in code ])

lbinary = compresso.raw_labels(cpso).tobytes()
print(len(lbinary))
print(len(binary))
# binary = lbinary + binary

with open("test-no-labels.bin", "wb") as f:
  f.write(binary)




