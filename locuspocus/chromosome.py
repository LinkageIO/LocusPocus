import reprlib
import numpy as np

from enum import Enum


class Nucleotide(Enum):
    A = 1
    a = 2
    C = 3
    c = 4
    G = 5
    g = 6
    T = 7
    t = 8
    U = 9
    u = 10
    N = 11
    n = 12


class Chromosome(object):
    """
    A Chromosome is a lightweight object which maps indices to
    string positions. It's a named sequence.
    NOTE: chromosome indices are different that python indices.
          Mainly, they are 1 indexed AND they are inclusive
          for start and stop positions.
    >>> from locuspocus import Chromosome
    >>> x = Chromosome('AAACCCTTTGGG')
    >>> x[1]
    'A'
    >>> x[1:5]
    'AAACC'
    >>> x[5:10]
    'CCTTTG'
    >>> len(x)
    12
    """

    def __init__(self, name, seq, *args):
        self.name = name
        # create the enumeration
        if isinstance(seq, np.ndarray):
            self.seq = seq
        else:
            self.seq = np.array([Nucleotide[x].value for x in seq])
        self._attrs = list(args)

    def __getitem__(self, pos):
        if isinstance(pos, slice):
            if pos.start < 1:
                raise ValueError("Genetic coordinates cannot start less than 1")
            return "".join(
                [Nucleotide(x).name for x in self.seq[max(0, pos.start - 1) : pos.stop]]
            )
        # chromosomes start at 1, python strings start at 0
        else:
            if pos < 1:
                raise ValueError("Genetic coordinates cannot start less than 1")
            return Nucleotide(self.seq[int(pos) - 1]).name

    def __len__(self):
        return len(self.seq)

    def __repr__(self):
        return "Chromosome({})".format(reprlib.repr("".join(self[1:100])))

    def __eq__(self, obj):
        if self.name == obj.name and all(self.seq == obj.seq):
            return True
        else:
            return False
