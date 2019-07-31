__version__ = '0.2.0'
__all__ = ['Locus','Loci']


#from .chromosome import Chromosome
#from .fasta import Fasta

#from .loci import Loci
#from .locus import Locus
from .loci import MemLoci as Loci
from .locus import MemLocus


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
    return MemLocus.create(chromosome,start,end,source,feature_type,strand,frame,name,attrs,parent,children)


