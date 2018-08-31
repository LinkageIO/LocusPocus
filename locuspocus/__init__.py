__version__ = '0.2.0'
__all__ = ['Locus','RefLoci','Fasta']

import pyximport
pyximport.install()

from .Chromosome import Chromosome
from .Fasta import Fasta

from .Locus import Locus
from .RefLoci import RefLoci

