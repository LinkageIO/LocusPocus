
import uuid

from . import LocusMixin
from locuspocus.locusdb import MemLocusDB

class MemLocus(LocusMixin,MemLocusDB):
    def __init__(self,UUID):
        # This super call just gets a handle to the in 
        # memory database storing all the current memloci
        super().__init__()
        self._UUID = UUID

    @classmethod
    def _from_UUID(cls,UUID):
        self = cls(UUID)
        return self

    @classmethod
    def create(
        cls,
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
        # Create a new UUID (convert to hex)
        self = cls(uuid.uuid4().hex)
        # Create a uuid
        # this starts a transaction
        with self._db:
            cur = self._db.cursor()
            # insert the core feature data
            cur.execute(
                '''
                INSERT INTO loci 
                    (UUID,chromosome,start,end,source,feature_type,strand,frame,name)
                    VALUES (?,?,?,?,?,?,?,?,?)
                ''',(
                # chrom/start/end are required, so cast them
                (self._UUID,str(chromosome),int(start),int(end),
                 source,feature_type,strand,frame,name))
            )

            # Add the key val pairs
            if attrs is not None:
                for key,val in attrs.items():
                    self[key] = val          

            # Handle Parent Child Relationships
            self.parent = parent
            self.children = children
            
            return self

    def __del__(self):
        self._db.cursor().execute('''
            DELETE FROM loci WHERE UUID = ?
        ''',(self._UUID,))
