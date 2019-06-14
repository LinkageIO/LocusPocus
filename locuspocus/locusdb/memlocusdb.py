import apsw as lite

from .locusdb import LocusDB

class MemLocusDB(LocusDB):

    _db = lite.Connection(':memory:')

    def __init__(self):
        # init tables defined in LocusDB class
        self._initialize_tables()

