
from locuspocus.locusdb import MemLocusDB
from .locimixin import LociMixin
from locuspocus.locus import MemLocus

class MemLoci(LociMixin,MemLocusDB):

    def __init__(self):
        super().__init__()

    def __getitem__(self, item):
        '''
            A convenience method to extract loci
            by their name attribute

            e.g.
            # Create a locus object
            >>> x = Locus('chr1',1,100,name='locus1')
            >>> y = Loci['locus1']
            >>> x == y
            True

        '''
        try:
            LID, = self._db.cursor().execute('''
                SELECT LID FROM loci WHERE name = ?
            ''',(item,)).fetchone()
        except TypeError:
            raise ValueError(f'{item} not in Loci')
        from locuspocus import Locus
        return memlocus(LID)


