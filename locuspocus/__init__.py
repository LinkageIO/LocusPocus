__version__ = '1.0.0-dev'
__all__ = ['Locus','Loci']


from .fasta import Fasta
from .chromosome import Chromosome
from .locus import MemLocus
from .loci import FrozenLoci


def Loci(name,basedir):
    return FrozenLoci(name, basedir)

def Locus(
        chromosome: str,
        start: int,
        end: int,

        source: str = 'locuspocus',
        feature_type: str = 'locus',
        strand: str = '+',
        frame: int = None,
        name: str = None,

        # Extra locus stuff
        attrs: dict = None,
        parent = None,
        children = None
):
    return MemLocus.create(
        chromosome,
        start,
        end,
        source,
        feature_type,
        strand,
        frame,
        name,
        attrs,
        parent,
        children
    )


