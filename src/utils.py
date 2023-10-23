# utils.py

def lookup(barcode: str, D: dict) -> int:
    """
    Lookup barcode in a dictionary.
    """
    s = D.get(barcode, 0) # return 0 for non-existent barcodes

    return s


def transpose(s, parity: int) -> int:
  """
  Transpose barcode given a parity.
  """
  t = (s-1)%parity+1

  return t


def encode(b: tuple[int], parity: int) -> int:
  """
  Encode recording based on barcode parity.
  """
  L = len(b)
  code = 0
  for i in range(L):
    s = b[i]
    t = transpose(s, parity)
    code += t*parity**(L-1-i)

  return code


def measure(code, parity: int) -> int:
  """
  Measure length of a recording.
  """
  Q = code
  L = 0
  while Q > 0:
    Q, R = divmod(Q, parity)
    if R == 0:
      Q -= 1
    L += 1

  return L