# bootstrap.py

"""Bootstrap a lineage tree from clonal DNA recordings."""

from src.utils import measure
import numpy as np

def vectorize(recordings: list) -> np.array:
  """
  Generate a probability vector for recordings within a clone.
  """
  D = {}
  for code in recordings:
    if code in D.keys():
      D[code] += 1
    else:
      D[code] = 1

  L = max(D.keys())+1
  v = np.zeros(L)
  for i in range(L):
    if i in D.keys():
      v[i] = D[i]

  v = np.array([x/v.sum() for x in v])

  return v


def generate_population_matrix(V: list) -> np.array:
  """
  Generate a population matrix from a list of clonal vectors.
  """
  L = max(v.size for v in V)
  P = np.zeros((len(V), L))
  for i, x in enumerate(V):
    P[i, :x.size] = x

  return P


def is_lineage_relation(i, j: int, c: np.array) -> bool:
  """
  Determine whether there is a lineage relation between two clones within a code vector.
  """
  for x, p in enumerate(c):
    if x == i or x == j:
      if p == 0:
        return False
    else:
      if p != 0:
        return False

  return True


def generate_lineage_matrix(clones: list, parity: int) -> np.array:
  """
  Generate a matrix of complexity-weighted distances of lineage recordings shared between clones within a population.
  """
  C = [vectorize(c) for c in clones]
  P = generate_population_matrix(C)
  L, E = P.shape
  A = np.zeros((L, L))
  for i in range(L):
    for j in range(L):
      if i == j:
        A[i][j] = 0
        continue
      S = 0
      for x in range(E):
        if is_lineage_relation(i, j, P[:,x]):
          s = 1-abs(P[i][x]-P[j][x])
          w = measure(x, parity)
          S += w*s
        else:
          continue
      A[i][j] = S
  
  return A

