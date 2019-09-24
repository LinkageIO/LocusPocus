import warnings

from locuspocus.locusdb import MemLocusDB
from .locimixin import LociMixin
from locuspocus.locus import MemLocus




class MemLoci(LociMixin,MemLocusDB):

    def __init__(self):
        warnings.warn(
            "Loci generated from MemLoci objects are not "
            "guaranteed to be stable as they are susceptible "
            "to undeterministic behavior due to the python "
            "garbage collector ",
            RuntimeWarning
        )
        super().__init__()

