__version__ = '1.0.0'
__all__ = ['Locus','RefLoci','Fasta']

import pyximport
pyximport.install()

from .Chromosome import Chromosome
from .Fasta import Fasta

from .Locus import Locus
from .RefLoci import RefLoci

from .Term import Term
