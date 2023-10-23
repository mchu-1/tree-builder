# profile.py

"""Profile clones based on DNA recordings."""

import yaml

def read_settings(config_f: str):
    """
    Generate settings dictionary from config YAML.
    """
    with open(config_f, "r") as f:
        print(f"Reading config file {config_f} ...")
        settings = yaml.safe_load(f)

    return settings


def get_sequences(input_f: str) -> list[str]:
    """
    Get sequences for a clone from a FASTQ file.
    """
    if input_f.split(".")[-1] != "fastq":
        raise ValueError("Input file should be in FASTQ format")

    sequences = []
    print(f"Accessing sequencing file {input_f} ...")
    with open(input_f, "r") as sequencing_f:
        marker = 0
        for line in sequencing_f:
            if marker == 1:
                new_sequence = line.rstrip()
                sequences.append(new_sequence)
                marker = 0
            elif line[0] == "@":
                marker = 1
                continue
            else:
                continue

    return sequences


def lookup(barcode: str, D: dict) -> int:
    """
    Lookup barcode in a dictionary.
    """
    s = D.get(barcode, 0) # return 0 for non-existent barcodes

    return s


def get_barcodes(sequence: str, spacer, h1, h2, s1, s2: str, l: int, D: dict) -> tuple[int]:
    """
    Get barcodes from a sequence.
    """
    b = []
    ix = sequence.find(spacer[:-3]) # index search at spacer stem
    if ix == -1:
        return tuple(b)
    else:
        ix += 17

    start_1, end_1 = h1[-l:], h2[:l]
    start_2, end_2 = h2[-l:], h1[:l]

    if sequence[ix: ix + 3] == s1:
        start, end = start_1, end_1
    elif sequence[ix: ix + 3] == s2:
        start, end = start_2, end_2
    else:
        return tuple(b)

    ix += 6
    while ix < len(sequence):
        if sequence[ix:][:l] != start:
            ix += 1
            continue
        elif sequence[ix + l + 4:][:l] != end:
            ix += 1
            continue
        else:
            new_barcode = sequence[ix + l:][:4]
            s = lookup(new_barcode, D)
            b.append(s)

            if start == start_1:
                start, end = start_2, end_2
            else:
                start, end = start_1, end_1

            ix += l

    b = tuple(b[::-1])

    return b


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


def get_recordings(input_f, config_f: str) -> list[int]:
    """
    Get recordings from a clone.
    """
    # Read sequences from input file
    sequences = get_sequences(input_f)

    # Read settings from config file
    settings = read_settings(config_f)
    spacer = settings["spacer"]
    h1, h2 = settings["h1"], settings["h2"]
    s1, s2 = settings["s1"], settings["s2"]
    l = settings["l"]
    D = settings["D"]
    parity = settings["parity"]

    # Get recordings
    recordings = []
    for sequence in sequences:
        b = get_barcodes(sequence, spacer, h1, h2, s1, s2, l, D)
        if len(b) == 0 or 0 in b: # remove zero-length recordings and unrecognized barcodes
            continue
        else:
            code = encode(b, parity)
            recordings.append(code)

    return recordings


