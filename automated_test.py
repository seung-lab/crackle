import crackle.crackcode
import numpy as np

def test_create_graph():
  labels = np.zeros((3,3), dtype=np.uint32)
  labels[1,1] = 1

  G = crackle.crackcode.create_graph(labels)

  assert set(G.nodes) == set([5,9,6,10])
  assert set(G.edges) == set([(5,9), (5,6), (6,10), (9,10)])

def test_create_crack_codes():
  labels = np.zeros((3,3), dtype=np.uint32)
  labels[1,1] = 1

  codes = crackle.crackcode.create_crack_codes(labels)
  ans = [[9, 0, 3, 1, 0, 2, 3, 3, 0]]
  assert codes == ans
