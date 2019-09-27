__version__ = '1.0.0'
__all__ = ['Locus','RefLoci','Fasta']

import pyximport
pyximport.install()

from .chromosome import Chromosome
from .fasta import Fasta

from .locus import Locus
from .loci import Loci

from .term import Term
