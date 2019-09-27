
from . import LocusMixin
from locuspocus.locusdb import FrozenLocusDB

class ThawedLocus(LocusMixin):
    def __init__(self,UUID,db):
        '''
            A thawed locus is created from a FrozenLoci object.
        '''
        super().__init__()
        self._UUID = UUID
        self._db = db
