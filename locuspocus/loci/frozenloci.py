from locuspocus.locusdb import FrozenLocusDB
from .locimixin import LociMixin
from locuspocus.locus import MemLocus

class FrozenLoci(LociMixin,FrozenLocusDB):

    def __init__(self,name,basedir=None):
        super().__init__(name,basedir=basedir)

