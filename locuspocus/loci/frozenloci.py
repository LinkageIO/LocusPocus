from locuspocus.locusdb import FrozenLocusDB
from .locimixin import LociMixin
from locuspocus.locus import ThawedLocus

class FrozenLoci(LociMixin,FrozenLocusDB):

    def __init__(self,name,basedir=None):
        super().__init__(name,basedir=basedir)

    def _get_locus_by_LID(self,LID):
        return ThawedLocus(LID,self._db)

    def add_locus():
        pass

    
