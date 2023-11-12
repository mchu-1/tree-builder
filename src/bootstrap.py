# bootstrap.py

"""Bootstrap a lineage tree from clonal DNA recordings."""

from src.utils import transit
import numpy as np

def generate_transtiion_matrix(recordings: list, parity: int) -> np.array:
  """
  Generate a transition matrix from a list of recordings within a clone.
  """
  L = max(len(r) for r in recordings)
  A = np.zeros((L+1, parity+2))
  for r in recordings:
    T = transit(r, parity)
    for t in T:
      A[t[0]][t[1]] += 1

  # Convert rows of A to probability distributions
  for i in range(L+1):
    A[i] = A[i]/sum(A[i])

  return A


def generate_lineage_matrix(clones: list, parity: int) -> np.array:
  """
  Generate lineage matrix by comparing transition matrices between clones within a population.
  """
  B = np.zeros(len(clones), len(clones))
  for i in range(len(clones)):
    for j in range(len(clones)):
      if i == j:
        B[i][j] = 0
        continue
      A_i = generate_transtiion_matrix(clones[i], parity)
      A_j = generate_transtiion_matrix(clones[j], parity)

      B[i][j] = np.linalg.norm(A_i-A_j) # Frobenius norm of transition matrices

  return B
