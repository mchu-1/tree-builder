# bootstrap.py

"""Bootstrap a lineage tree from clonal DNA recordings."""

from utils import transit
import numpy as np

def generate_transition_matrix(recordings: list, parity: int) -> np.array:
  """
  Generate a transition matrix from a list of recordings within a clone.
  """
  L = max(len(r) for r in recordings)
  A = np.zeros((L+1, (parity+2)**2))
  for r in recordings:
    T = transit(r, parity)
    a = np.zeros((parity+2, parity+2))
    for i, t in enumerate(T):
      A[i][t[0]*(parity+2)+t[1]] += 1

  # Convert rows of A to probability distributions
  for i in range(L+1):
    A[i] = A[i]/sum(A[i])

  return A


def generate_population_matrix(clones: list, parity: int) -> np.array:
  """
  Generate population matrix of transitions for all clones.
  """
  # Generate transition matrices for all clones
  A = [generate_transition_matrix(clone, parity) for clone in clones]
  # Flatten transition matrices
  A = [a.flatten() for a in A]
  dim = max(len(a) for a in A)
  # Pad all transition matrices to the same dimension with zeros
  A = [np.pad(a, (0, dim-len(a)), 'constant') for a in A]
      
  return A


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
  Generate a matrix of complexity-weighted similarity between clonal transitions within a population.
  """
  P = generate_population_matrix(clones, parity)
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
          w = x//(parity+2)
          S += (parity**w)*s
        else:
          continue
      A[i][j] = S
  
  return A
