import apsw as lite

from .locusdbmixin import LocusDBMixin

class MemLocusDB(LocusDBMixin):

    _db = lite.Connection(':memory:')

    def __init__(self):
        # init tables defined in LocusDB class
        self._initialize_tables()

