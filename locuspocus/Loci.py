#!/usr/bin/python3
import itertools
import random
import gzip
import logging

import numpy as np
import scipy as sp
import apsw as lite

from minus80 import Freezable
from collections.abc import Iterable
from functools import lru_cache,wraps
from contextlib import contextmanager
from typing import List

from locuspocus import Locus
from locuspocus.Exceptions import ZeroWindowError,MissingLocusError,StrandError

__all__ = ['Loci']

# --------------------------------------------------
#
#       Decorators
#
# --------------------------------------------------

#def invalidates_primary_loci_cache(fn):
#    '''
#    Decorate any methods that invalidate the list
#    of primary loci. (e.g. the _primary_LIDS method)
#    '''
#    @wraps(fn)
#    def wrapped(self,*args,**kwargs):
#        retval = fn(self,*args,**kwargs) 
#        self._primary_LIDS.cache_clear()
#        return retval
#    return wrapped

def accepts_loci(fn):
    '''
    This decorator augments methods that take as their first
    argument a Locus object. It allows the method to also accept
    an iterable of Locus objects and maps the method to the
    Locus objects in the iterable.
    '''
    @wraps(fn)
    def wrapped(self,loci,*args,**kwargs):
        if not isinstance(loci,Locus):
            return [fn(self,l,*args,**kwargs) for l in loci]
        else:
            return fn(self,loci,*args,**kwargs)
    return wrapped

class Loci(Freezable):
    '''
        Loci are more than the sum of their parts. They have a name and
        represent something bigger than theirselves. They are important. They
        live on the disk in a database.
    '''
    _db = lite.Connection(':memory:')
    # Set up a logger
    log = logging.getLogger(__name__)
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(asctime)s %(name)-12s %(levelname)-8s %(message)s")
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(logging.INFO)

    def __init__(self, name=None, basedir=None):
        self.name = name
        if self.name is not None:
            # this sets up the Freezable API which overloads the
            # class attribute: self._db
            super().__init__(name, basedir=basedir)
        # initialize some instance variables
        self._initialize_tables()
        self.primary_type = 'locus'

    def add_locus(
        self,

        chromosome: str,
        start: int,
        end: int,

        source: str,
        feature_type: str,
        strand: str,
        frame: int,
        name: str,

        attrs: dict = None,
        parent: Locus = None,
        children: List = None,

        cur = None
    ) -> int:
        '''
            Add locus to the database. 

            Parameters
            ----------
            cur : a db cursor
                An optional cursor object to use. If none, a 
                cursor will be created. If importing many 
                loci in a loop, use a bulk transaction.

            Returns
            -------
            The locus ID (LID) of the freshly added locus
        '''
        if cur is None:
            cur = self._db.cursor()
        # this starts a transaction
        with self._db:
            # insert the core feature data
            cur.execute(
                '''
                INSERT INTO loci 
                    (chromosome,start,end,source,feature_type,strand,frame,name)
                    VALUES (?,?,?,?,?,?,?,?)
                ''',
                (chromosome,start,end,source,feature_type,strand,frame,name)
            )
            # get the fresh LID
            (LID,) = cur.execute('SELECT last_insert_rowid()').fetchone()
            if LID is None: #pragma: no cover
                # I dont know when this would happen without another exception being thrown
                raise ValueError(f"{locus} was not assigned a valid LID!")
            # Add the key val pairs
            if attrs is not None:
                for key,val in attrs.items():
                    cur.execute(
                        '''
                        INSERT INTO loci_attrs 
                        (LID,key,val) 
                        VALUES (?,?,?)
                        ''', 
                        (LID,key,val)
                    )
             # Check to see if the locus is a primary feature
            if feature_type == self.primary_type:
                cur.execute('''
                    INSERT INTO primary_loci (LID) VALUES (?)
                ''',(LID,))

            # Handle Parent Child Relationships
            if parent is not None:
                try:
                    cur.execute('''
                        INSERT INTO relationships (parent,child)
                        VALUES (?,?)
                    ''',(parent._LID,LID))
                except AttributeError: 
                    raise ValueError('Parent of locus must be a locus')
            if children is not None:
                try:
                    children_LIDs = [c._LID for c in locus.subloci]
                except AttributeError as e:
                    raise ValueError('All children of locus must also be of type Locus')
                # Add relationships  
                for CID in children_LIDs:
                    # children can only have one parent
                    # detach any children and reassign as parent
                    cur.execute('''
                        DELETE FROM relationships
                        WHERE child IN
                    ''')
                    cur.execute(
                        '''
                        INSERT INTO relationships (parent,child)
                        VALUES (?,?)
                        ''',
                        (LID,CID)
                    )
            # Add the position to the R*Tree
            cur.execute(
                '''
                INSERT INTO positions (LID,start,end,chromosome) VALUES (?,?,?,?)
                ''',
                (LID,start,end,chromosome)
            )

        return LID
    
    def del_locus(self,locus: Locus) -> bool:
        '''
            Remove a locus from the set of Loci
        '''
        if locus._ref != self:
            raise ValueError('Locus not in database!')
        self._db.cursor().execute('''
            DELETE FROM loci WHERE LID = ?
        ''',(locus._LID,))

    @property
    def _primary_LIDS(self) -> List[int]:
        '''
        A cached list of primary Locus IDs (LIDS) available in the Loci database.
        This list can change if the primary locus changes, for example through
        the `Loci.set_primary_feature_type()` method. 
        
        Note: this method takes no arguments.
        '''
        LIDS = [x[0] for x in self._db.cursor().execute('SELECT LID FROM primary_loci')]
        return LIDS

    def __len__(self) -> int:
        '''
        Returns the number of primary loci in the dataset. Gets called by the len()
        built-in function

        >>> len(ref)
        42
        '''
        l, = self._db.cursor().execute(
            'SELECT COUNT(*) FROM primary_loci'
        ).fetchone()
        return l

    def _get_locus_by_LID(self,LID: int) -> Locus:
        '''
        Get a locus by its LID
        
        Parameters
        ----------
        LID : int
            A Locus ID. These are assigned to Locus objects when
            they are added to the Loci database.

        Returns
        -------
        The Locus.

        Raises
        ------
        `MissingLocusError` if there is no Locus in the database with that LID.
        '''
        lid_exists, = self._db.cursor().execute('SELECT COUNT(*) FROM loci WHERE LID = ? ',(LID,)).fetchone()
        if lid_exists == 0:
            raise MissingLocusError
        return Locus(LID,self)


    def __getitem__(self, item):
        '''
            A convenience method to extract loci from the reference genome.
        '''
        raise NotImplementedError

    def __iter__(self):
        return (self._get_locus_by_LID(l) for l in self._primary_LIDS())

    @contextmanager
    def filter_feature_type(self,feature_type):
        '''
        Switch primary feature types within a context block.
        Useful for temporary switching in a with block.
        
        >>> ref = locuspocus.Loci('test')
        >>> with ref.filter_feature_type('exon'):
                ...
                ... 'Do something with exons here'
                ...
        >>> 'primary type is now back to whatever it was before' 
        '''
        cur_lids = self._primary_LIDS()
        self.set_primary_feature_type(feature_type,clear_previous=True)
        yield
        with self._bulk_transaction() as cur:
            cur.execute('DELETE FROM primary_loci')
            cur.executemany('''
                INSERT INTO primary_loci (LID) 
                VALUES (?)
                ''', ((x,) for x in cur_lids)
            )
        # reset the LID cache here
        self._primary_LIDS.cache_clear()


    def set_primary_feature_type(self,feature_type,clear_previous=True):
        '''
        Set the primary feature type.
        
        Parameters
        ----------
        feature_type : str
            The feature type of the primary locus type
        clear_previous : bool (default: True)
            If true, delete previous entries for primary loci.
            If false, append types

        '''
        with self._bulk_transaction() as cur: 
            if clear_previous:
                cur.execute('DELETE FROM primary_loci')
            cur.execute('''
                INSERT INTO primary_loci (LID) 
                SELECT LID FROM loci WHERE feature_type = ?
            ''',(feature_type,))


    def rand(self, n=1, distinct=True, autopop=True):
        '''
            Fetch random Loci

            Parameters
            ----------
            n : int (default=1)
                The number of random locus objects to fetch
            distinct : bool (default=True)
                If True, the loci are all guaranteed to be distinct, 
                AKA, sampling without replacement. If False, the
                sampling is with replacement so there may be
                duplucates in the out output, especially if the number
                of n is high.
            autopop : bool (default: True)
                If true and only 1 locus is requested, a Locus object
                will be returned instead of a list (with a single element)

            Returns
            -------
            A list of n Locus objects

        '''
        import random
        loci = set()
        LIDS = self._primary_LIDS()
        if n > len(LIDS):
            raise ValueError('More than the maximum loci in the database was requested')
        if distinct == True:
            LIDs = random.sample(self._primary_LIDS(),n)
        else:
            LIDS = random.choices(self._primary_LIDS(),n)
        loci = [self._get_locus_by_LID(x) for x in LIDs]
        if autopop and len(loci) == 1:
            loci = loci[0]
        return loci

    def feature_types(self, print_tree=True): #pragma: no cover
        '''
        Returns a summary of the feature types represented in 
        the Loci database
        
        Parameters
        ----------
        print_tree : bool (default: True)
            If True, prints the result before returning. 

        Returns
        -------
        An anytree Node object containing the root node.

        '''
        raise NotImplementedError('This method is BUGGY')
        from anytree import Node, RenderTree

    @accepts_loci 
    def within(
        self, 
        locus, 
        partial=False, 
        ignore_strand=False,
        same_strand=False
    ):
        '''
        Returns the Loci that are within the start/stop boundaries
        of an input locus.

        By default, the entire locus needs to be within the
        coordinates of the input locus, however, this can be
        toggled with the `partial` argument.
        
        NOTE: this loci are returned in order of 3' to 5'
              based on the strand of the input locus. This
              bahavior can be changed by using the `ignore_strand` 
              option.

        __________________Ascii Example___________________________
        These features get returned (y: yes, n:no)

        partial=False (default):
            nnnnn yyyy yyy   nnnnnnnnnnn    
        partial=True:
            yyyyy yyyy yyy   yyyyyyyyyyy    
        __________________________________________________________

              start           end
        -------[===============]-------------

        The loci will be returned based on the strand of the input
        locus. The following numbers represent the indices in the 
        generator returned by the method.

        (+) stranded Locus (assume partial=True):
              0    1    2       3
            yyyyy yyyy yyy   yyyyyyyyyyy    

        (-) stranded Locus (assume partial=True):
              3    2    1       0
            yyyyy yyyy yyy   yyyyyyyyyyy    
        __________________________________________________________

        Parameters
        ----------
        locus : Locus
            A locus defining the genomic interval to 
            extract loci within.
        partial : bool (default=False)
            When True, include loci that partially overlap.
            See example ASCII above.
        ignore_strand : bool (default=False)
            If ignore_strand is True, the method will
            assume a (+) strand, otherwise it will 
            respect the strand of the input locus
            in relation to the order the resulting
            Loci are returned. See Ascii example.
            Note: this cannot be True if same_strand
            is also set to True, an exception will be 
            raised.
        same_strand : bool (default: False)
            If True, only Loci on the same strand
            as the input locus will be returned,
            otherwise, the method will return loci
            on either strand. Note: this cannot be 
            set if `ignore_strand` is also true. An
            exception will be raised.
        '''
        if ignore_strand and same_strand:
            raise ValueError('`ignore_strand` and `same_strand` cannot both be True')
        # set up variables to use based on 'partial' flag
        cur = self._db.cursor()
        # Calculate the correct strand orientation
        if locus.strand == '+' or ignore_strand == True:
            if partial == False:
                anchor = f'l.start > {locus.start} AND l.end < {locus.end}'
                order = 'l.start ASC'
                index = 'locus_start'
            else:
                anchor = f'l.end > {locus.start} AND l.start < {locus.end}'
                order = 'l.end ASC'
                index = 'locus_end'
        elif locus.strand == '-':
            if partial == False:
                anchor = f'l.end < {locus.end} AND l.start > {locus.start}'
                order = 'l.end DESC'
                index = 'locus_end'
            else:
                anchor = f'l.start < {locus.end} AND l.end  > {locus.start}'
                order = 'l.start DESC'
                index = 'locus_start'
        else:
            raise StrandError
        query = f'''
            SELECT l.LID FROM primary_loci p, loci l 
            INDEXED BY {index} 
            WHERE l.chromosome = '{locus.chromosome}' 
            AND {anchor} 
            AND l.LID = p.LID 
            ORDER BY {order}; 
        '''
        LIDS = cur.execute(query)
        for x, in LIDS:
            l = self._get_locus_by_LID(x) 
            if same_strand == True and l.strand != locus.strand:
                continue
            yield l

    @accepts_loci
    def upstream_loci(
        self, 
        locus, 
        n=np.inf, 
        max_distance=10e100, 
        partial=False,
        same_strand=False,
        force_strand=None
    ):
        '''
            Find loci upstream of a locus.

            Loci are ordered so that the nearest loci are
            at the beginning of the list.

            NOTE: this respects the strand of the locus unless
                  `force_strand` is set.

            __________________Ascii Example___________________________
            These features get returned (y: yes, n:no)

            partial=False
            nnn  nnnnnnn   nnnnnn   yyyyy  nnnnnn nnnn nnnn nnnnnnnn

            partial=True
            nnn  nnnnnnn   yyyyyy   yyyyy  yyyyyy nnnn nnnn nnnnnnnn
            __________________________________________________________
                                            
            (+) stranded Locus:            start             end
                -----------------------------[================]--
                             <_______________| upstream

            (-) stranded Locus:         
            start          end
             [===============]----------------------------------
                   upstream  |_______________> 

            __________________________________________________________

            Parameters
            ----------
                locus : Locus object
                    The locus object for which to fetch the upstream
                    loci.
                n : int (default: infinite)
                    This maximum number of Loci to return. Note that
                    there are cases where fewer loci will be available,
                    e.g.: at the boundaries of chromosomes, in which
                    case fewer loci will be returned.
                max_distance : float (deflocus.stranded_startp.inf)
                    The maximum distance 
                partial : bool (default: False)
                    A flag indicating whether or not to return
                    loci that are partially within the boundaries
                    of the input Locus. See Ascii example above.
                same_strand : bool (default: False)
                    If True, only Loci on the same strand
                    as the input locus will be returned,
                    otherwise, the method will return loci
                    on either strand.
        '''
        # calculate the start and stop anchors 
        start,end = sorted([locus.stranded_start, locus.upstream(max_distance)])
        # The dummy locus needs to have the opposite "strand" so the loci
        # are returned in the correct order
        dummy_strand = '+' if locus.strand == '-' else '-'
        # create a dummy locus
        upstream_region = Locus(locus.chromosome,start,end,strand=dummy_strand)
        # return loci within dummy locus coordinates
        loci = self.within(
            upstream_region, 
            partial=partial
        )
        # Keep track of how many values have been yielded 
        i = 1
        for x in loci:
            # If we've yielded enough values, stop
            if i > n:
                break
            # Check to see if if we should filter loci
            if same_strand and x.strand != locus.strand:
                continue
            yield x
            i+=1

    @accepts_loci
    def downstream_loci(
        self, 
        locus, 
        n=np.inf, 
        max_distance=10e100, 
        partial=False,
        ignore_strand=False,
        same_strand=False
    ):
        '''
            Returns loci downstream of a locus. 
            
            Loci are ordered so that the nearest loci are 
            at the beginning of the list.

            NOTE: this respects the strand of the locus unless
                  `ignore_strand` is set to True
            
            __________________Ascii Example___________________________
            These features get returned (y: yes, n:no)

            partial=False
            nnn  nnnnnnn   nnnnnn   yyyyy  nnnnnn nnnn nnnn nnnnnnnn

            partial=True
            nnn  nnnnnnn   yyyyyy   yyyyy  yyyyyy nnnn nnnn nnnnnnnn
            __________________________________________________________
                                            
            (+) stranded Locus:         
            start          end
             [===============]----------------------------------
                 downstream  |_______________> 

            (-) stranded Locus:            start             end
                -----------------------------[================]--
                             <_______________| downstream

            __________________________________________________________



            Parameters
            ----------
                locus : Locus object
                    The locus object for which to fetch the upstream
                    loci.
                n : int (default: infinite)
                    This maximum number of Loci to return. Note that
                    there are cases where fewer loci will be available,
                    e.g.: at the boundaries of chromosomes, in which
                    case fewer loci will be returned.
                max_distance : float (deflocus.stranded_startp.inf)
                    The maximum distance 
                partial : bool (default: False)
                    A flag indicating whether or not to return
                    loci that are partially within the boundaries
                    of the input Locus. See Ascii example above.
                same_strand : bool (default: False)
                    If True, only Loci on the same strand
                    as the input locus will be returned,
                    otherwise, the method will return loci
                    on either strand.
        '''
        # calculate the start and stop anchors 
        start,end = sorted([locus.stranded_end, locus.downstream(max_distance)])
        # The dummy locus needs to have the same strand so that
        # are returned in the correct order
        dummy_strand = '+' if locus.strand == '+' else '-'
        # create a dummy locus
        downstream_region = Locus(locus.chromosome,start,end,strand=dummy_strand)
        # return loci within dummy locus coordinates
        loci = self.within(
            downstream_region, 
            partial=partial
        )
        # Keep track of how many values have been yielded 
        i = 1
        for x in loci:
            # If we've yielded enough values, stop
            if i > n:
                break
            # Check to see if if we should filter loci
            if same_strand and x.strand != locus.strand:
                continue
            yield x
            i+=1

    def flanking_loci(
        self,
        locus,
        n=np.inf,
        max_distance=10e100,
        partial=False,
        same_strand=False
    ):
        '''
            This is a convenience method to return loci both upstream 
            and downstream of an input locus. The options and flags are
            simply passed to the aforementioned function.

            Parameters
            ----------
            locus : Locus object
                The input locus for which to return flanking loci
            n : int (default=infinite)
                
        '''
        kwargs = {
            'n':n,
            'max_distance':max_distance,
            'partial':partial,
            'same_strand':same_strand
        }
        return (
            self.upstream_loci(
                locus,
                **kwargs
            ),
            self.downstream_loci(
                locus,
                **kwargs
            )
        )

    def encompassing_loci(self, locus):
        '''
            Returns the Loci encompassing the locus. In other words
            if a locus (e.g. a SNP) is inside of another locus, i.e. the
            start of the locus is upstream and the end of the locus
            is downstream of the locus boundaries, this method will
            return it. 

            Parameters
            ----------
            locus : Locus object

            Returns
            -------
            Loci that encompass the input loci
        '''
        cur = self._db.cursor()
        LIDS = cur.execute('''
            SELECT l.LId FROM positions p, primary_loci l
            WHERE p.chromosome = ?
            AND p.start < ?
            AND p.end > ?
            AND p.LID = l.LID
        ''',(locus.chromosome,locus.start,locus.end))
        for x, in LIDS:
            yield self._get_locus_by_LID(x) 

    def _nuke_tables(self):
        cur = self._db.cursor()
        cur.execute(
            '''
            DROP TABLE IF EXISTS loci;
            DROP TABLE IF EXISTS loci_attrs;
            DROP TABLE IF EXISTS aliases;
            DROP VIEW IF EXISTS named_loci;
            DROP TABLE IF EXISTS relationships;
            DROP TABLE IF EXISTS positions;
            DROP TABLE IF EXISTS primary_loci;
            '''
        )
        self._initialize_tables()

    def import_gff(
        self, 
        filename: str, 
        feature_type="*", 
        ID_attr="ID", 
        parent_attr='Parent',
        attr_split="="
    ) -> None:
        '''
            Imports Loci from a gff (General Feature Format) file.
            See more about the format here:
            http://www.ensembl.org/info/website/upload/gff.html

            Parameters
            ----------

            filename : str
                The path to the GFF file.
            ID_attr : str (default: ID)
                The key in the attribute column which designates the ID or
                name of the feature.
            parent_attr : str (default: Parent)
                The key in the attribute column which designates the Parent of
                the Locus
            attr_split : str (default: '=')
                The delimiter for keys and values in the attribute column
        '''
        self.log.info(f"Importing Loci from {filename}")
        if filename.endswith(".gz"):
            IN = gzip.open(filename, "rt")
        else:
            IN = open(filename, "r")
        loci = []
        name_index = {}
        total_loci = 0
        for i,line in enumerate(IN):
            total_loci += 1
            # skip comment lines
            if line.startswith("#"):
                continue
            # Get the main information
            (
                chromosome,
                source,
                feature,
                start,
                end,
                score,
                strand,
                frame,
                attributes,
            ) = line.strip().split("\t")
            # Cast data into appropriate types
            start = int(start)
            end = int(end)
            strand = None if strand == '.' else strand
            frame = None if frame == '.' else int(frame)
            # Get the attributes
            attributes = dict(
                [
                    (field.strip().split(attr_split))
                    for field in attributes.strip(";").split(";")
                ]
            )
            # Store the score in the attrs if it exists
            if score != '.':
                attributes['score'] = float(score)
            # Parse out the Identifier
            if ID_attr in attributes:
                name = attributes[ID_attr]
                #del attributes[ID_attr]
            else:
                name = None
            # Parse out the parent info
            if parent_attr in attributes:
                parent = attributes[parent_attr]
                #del attributes[parent_attr]
            else:
                parent = None
            l =  Locus(
                    chromosome=chromosome, 
                    start=start, 
                    source=source, 
                    feature_type=feature, 
                    end=end, 
                    strand=strand, 
                    frame=frame, 
                    name=name, 
                    attrs=attributes, 
                )
            if name is not None:    
                name_index[name] = l
            if parent is None:
                loci.append(l)
            else:
                name_index[parent].add_sublocus(l)
        self.log.info((
            f'Found {len(loci)} top level loci and {total_loci} '
            f'total loci, adding them to database'
        ))
        IN.close()
        with self._bulk_transaction() as cur:
            for l in loci:
                self.add_locus(l,cur=cur)
        self.log.info('Done!')
        return None


    def _initialize_tables(self):
        '''
        Initializes the Tables holding all the information
        about the Loci.
        '''
        cur = self._db.cursor()
        cur.execute('''
            PRAGMA foreign_keys = ON;
        ''')
        cur.execute(
            '''
            CREATE TABLE IF NOT EXISTS loci (
                LID INTEGER PRIMARY KEY AUTOINCREMENT,
                
                /* Store the locus values  */
                chromosome TEXT NOT NULL,
                start INTEGER NOT NULL,
                end INTEGER,

                source TEXT,
                feature_type TEXT,
                strand TEXT,
                frame INT,

                /* Store things to make my life easier */
                name TEXT
                
            );
            '''
        )

        cur.execute('''
            CREATE INDEX IF NOT EXISTS locus_id ON loci (name);
            CREATE INDEX IF NOT EXISTS locus_chromosome ON loci (chromosome);
            CREATE INDEX IF NOT EXISTS locus_start ON loci (start);
            CREATE INDEX IF NOT EXISTS locus_feature_type ON loci (feature_type);
            CREATE INDEX IF NOT EXISTS locus_end ON loci (end);
            '''
        )
        cur.execute(
        # Create a table that contains loci attribute mapping
        '''
            CREATE TABLE IF NOT EXISTS loci_attrs (
                LID INT REFERENCES loci(LID) ON DELETE CASCADE,
                key TEXT,
                val TEXT,
                UNIQUE(LID,key)
            );
            '''
        )
        cur.execute(
        # Create a table with parent-child relationships        
        '''
            CREATE TABLE IF NOT EXISTS relationships (
                parent INT REFERENCES loci(LID) ON DELETE CASCADE,
                child INT REFERENCES loci(LID) ON DELETE CASCADE
            );
            CREATE INDEX IF NOT EXISTS relationships_parent ON relationships (parent);
            CREATE INDEX IF NOT EXISTS relationships_child ON relationships (child);
        '''
        )
        # Create a R*Tree table so we can efficiently query by ranges
        cur.execute(
        '''
            CREATE VIRTUAL TABLE IF NOT EXISTS positions USING rtree_i32( 
                LID, 
                start INT,
                end INT,
                +chromosome TEXT
            );
        '''
        )
        cur.execute(
        '''
            CREATE TABLE IF NOT EXISTS primary_loci (
                LID INT PRIMARY KEY REFERENCES loci(LID) ON DELETE CASCADE
            );
        '''
        )



