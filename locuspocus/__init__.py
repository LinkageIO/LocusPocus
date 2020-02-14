__version__ = '1.0.2'

__all__ = ['Locus','RefLoci','Fasta']

import logging

from .chromosome import Chromosome
from .fasta import Fasta

from .locus import Locus
from .loci import Loci

from .term import Term

log = logging.getLogger('locuspocus')
# One of DEBUG, INFO, WARNING, ERROR, CRITICAL
log.setLevel(logging.INFO)
# Set up the console handler
ch = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s - %(message)s',datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
ch.setLevel(logging.INFO)
log.addHandler(ch)


