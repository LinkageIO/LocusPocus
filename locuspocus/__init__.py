__version__ = '1.0.0-dev'
__all__ = ['Locus','Loci']


from .locus import MemLocus
from .loci import Loci



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


