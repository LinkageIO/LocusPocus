from minus80 import Freezable

from .locusdbmixin import LocusDBMixin

class FrozenLocusDB(LocusDBMixin, Freezable):

    def __init__(self,name,basedir=None):
        super().__init__(name,basedir=basedir)
        # Ititialize the minus80 API
        self._db = self.m80.db
        self._initialize_tables()
