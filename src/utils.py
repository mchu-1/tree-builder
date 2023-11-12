# utils.py

def lookup(barcode: str, D: dict) -> int:
    """
    Lookup barcode in a dictionary.
    """
    s = D.get(barcode, -1) # return -1 for unknown barcodes
    return s


def transpose(s, parity: int) -> int:
  """
  Transpose barcode given a parity.
  """
  t = (s-1)%parity+1

  return t


def transit(b: tuple[int], parity: int) -> list:
  """
  Get all barcode transitions in a recording.
  """
  L = len(b)
  T = [0] # placeholder for beginning of recording
  for i in range(L):
    s = b[i]
    t = transpose(s, parity)
    T.append(t)
  T.append(parity+1) # placeholder for end of recording

  T = list(zip(T[:-1], T[1:]))

  return T