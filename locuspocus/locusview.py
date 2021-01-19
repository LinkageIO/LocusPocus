
from .locus import Locus,LocusAttrs,SubLoci

__all__ = ['LocusView']

class AttrsView(LocusAttrs):
    def __init__(self, parent):
        self.parent = parent

    @property
    def empty(self):
        if len(self) == 0:
            return True
        return False

    @property
    def table(self):
        if self.parent.is_sublocus:
            return 'subloci_attrs'
        else:
            return 'loci_attrs'

    def __len__(self):
        cur = self.parent._ref.m80.db.cursor()
        cur.execute(f'''
            SELECT COUNT(*) FROM {self.table}
            WHERE LID = ?
            ''',
            (self.parent._LID,)
        )
        return cur.fetchone()[0]

    def keys(self):
        cur = self.parent._ref.m80.db.cursor()
        results = cur.execute(f'''
            SELECT key FROM {self.table}
            WHERE LID = ?
            ''',
            (self.parent._LID,)
        )
        return [k[0] for k in results]

    def values(self):
        cur = self.parent._ref.m80.db.cursor()
        results = cur.execute(f'''
            SELECT val FROM {self.table} 
            WHERE LID = ?
            ''',
            (self.parent._LID,)
        ) 
        return [k[0] for k in results]

    def items(self):
        return zip(self.keys(),self.values())

    def __getitem__(self,key):
        cur = self.parent._ref.m80.db.cursor()
        try:
            val, = cur.execute(f'''
                SELECT val FROM {self.table}
                WHERE LID = ? AND key = ?
            ''',(self.parent._LID,key)).fetchone()
        except TypeError:
            raise KeyError(f'"{key}" in in attrs')
        return val 

    def __setitem__(self,key,val):
        cur = self.parent._ref.m80.db.cursor()
        cur.execute(f'''
            INSERT OR REPLACE INTO {self.table}
            (LID,key,val)
            VALUES (?,?,?)
        ''',(self.parent._LID,key,val))

    def __repr__(self):
        return '{'+','.join([':'.join([x,y]) for x,y in self.items()]) +'}'


class SubLociView(SubLoci):
    # A restricted list interface to subloci
    def __init__(self, parent):
        self.parent = parent 

    @property
    def empty(self):
        if len(self) == 0:
            return True
        return False
   
    @property
    def _LID_query(self):
        if not self.parent.is_sublocus:
            query = '''
                SELECT LID FROM subloci 
                WHERE root_LID = ? 
                AND parent_LID IS NULL
            '''
        else:
            query = '''
                SELECT LID FROM subloci 
                WHERE parent_LID = ?
            '''
        return query

    def __iter__(self):
        cur = self.parent._ref.m80.db.cursor()
        LIDs = (LID for (LID,) in cur.execute(
            self._LID_query,
            (self.parent._LID,)
        ))
        return (LocusView(x, self.parent._ref, sublocus=True) for x in LIDs)

    def add(self,locus):
        raise NotImplementedError

    def __getitem__(self,index):
        query = self._LID_query + 'ORDER BY LID LIMIT 1 OFFSET ?'
        (LID,) = self.parent._ref.m80.db.cursor().execute(
            query,
            (self.parent._LID,index)
        ).fetchone()
        return LocusView(
            LID,
            self.parent._ref,
            sublocus=True
        )

    def __len__(self):
        return self.parent._ref.m80.db.cursor().execute(
            self._LID_query.replace(
                'SELECT LID FROM subloci',
                'SELECT COUNT(LID) FROM subloci'
            ),
            (self.parent._LID,)
        ).fetchone()[0]

    def __repr__(self):
        if self.empty:
            return '[]'
        return '[\n'+'\t\n'.join([repr(x) for x in self])+'\n]'

class LocusView(Locus):
    '''
        A LocusView is an efficient way to acces
        Locus objects stored in a Loci database 
    '''
    def __init__(
        self, 
        LID: int, 
        refloci: "Loci", 
        sublocus: bool = False
    ):
        self._LID = LID
        self._ref = refloci
        self._sublocus = sublocus
        self.attrs = AttrsView(self)
        self.subloci = SubLociView(self) 

    @property
    def is_sublocus(self):
        if self._sublocus == True:
            return True
        else:
            return False

    @property
    def table(self):
        if self.is_sublocus:
            return 'subloci'
        else:
            return 'loci'

    def _property(self,name):
        val, = self._ref.m80.db.cursor().execute(
            f"SELECT {name} FROM {self.table} WHERE LID = ?",
            (self._LID,)
        ).fetchone()
        return val

    @property
    def chromosome(self):
        return self._property('chromosome')
        
    @property
    def start(self):
        return self._property('start')

    @property
    def end(self):
        return self._property('end')

    @property
    def feature_type(self):
        return self._property('feature_type')

    @property
    def strand(self):
        return self._property('strand')

    @property
    def frame(self):
        return self._property('frame')

    @property
    def name(self):
        return self._property('name')

    @property
    def source(self):
        return self._property('source')

    def add_sublocus(self,locus):
        raise NotImplementedError
