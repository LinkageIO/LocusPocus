#!/usr/bin/python3
import itertools
import random
import scipy as sp
import gzip
import logging

import numpy as np

from minus80 import Freezable
from collections.abc import Iterable
from functools import lru_cache,wraps
from contextlib import contextmanager
from typing import List

from .locus import Locus
from .locusview import LocusView
from .exceptions import ZeroWindowError,MissingLocusError,StrandError

__all__ = ['RefLoci']


# --------------------------------------------------
#
#       Decorators
#
# --------------------------------------------------

def invalidates_primary_loci_cache(fn):
    '''
    Decorate any methods that invalidate the list
    of primary loci. (e.g. the _primary_LIDS method)
    '''
    @wraps(fn)
    def wrapped(self,*args,**kwargs):
        retval = fn(self,*args,**kwargs) 
        self._primary_LIDS.cache_clear()
        return retval
    return wrapped

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
        RefLoci are more than the sum of their parts. They have a name and
        represent something bigger than theirselves. They are important. They
        live on the disk in a database.
    '''

    def __init__(self, name, basedir=None):
        # set up the freezable API
        super().__init__(name, basedir=basedir)
        self.name = name
        self._initialize_tables()
        # Set up a logger
        self.log = logging.getLogger(__name__)
        handler = logging.StreamHandler()
        formatter = logging.Formatter("%(asctime)s %(name)-12s %(levelname)-8s %(message)s")
        handler.setFormatter(formatter)
        self.log.addHandler(handler)
        self.log.setLevel(logging.INFO)

    @lru_cache(maxsize=1)
    def _primary_LIDS(self) -> List[int]:
        '''
        A cached list of primary Locus IDs (LIDS) available in the RefLoci database.
        This list can change if the primary locus changes, for example through
        the `RefLoci.set_primary_feature_type()` method. 
        
        Note: this method takes no arguments.
        '''
        self.log.info("Caching LIDS from database")
        LIDS = [x[0] for x in self.m80.db.cursor().execute('SELECT LID FROM primary_loci')]
        return LIDS


    def __len__(self) -> int:
        '''
        Returns the number of primary loci in the dataset. Gets called by the len()
        built-in function

        >>> len(ref)
        42
        '''
        return len(self._primary_LIDS())

    def _get_locus_by_LID(self,LID: int) -> LocusView:
        '''
        Get a locus by its LID
        
        Parameters
        ----------
        LID : int
            A Locus ID. These are assigned to Locus objects when
            they are added to the RefLoci database.

        Returns
        -------
        The LocusView.

        Raises
        ------
        `MissingLocusError` if there is no Locus in the database with that LID.
        '''
        lid_exists, = self.m80.db.cursor().execute('SELECT COUNT(*) FROM loci WHERE LID = ? ',(LID,)).fetchone()
        if lid_exists == 0:
            raise MissingLocusError
        return LocusView(LID,self)

    def _get_LID(self,locus: Locus) -> int: #pragma: no cover
        '''
            Return the Locus Identifier used internally by RefLoci

            Parameters
            ----------
            locus : one of (str,Locus)
                The locus for which to find the LID of, this can
                be either a Locus object OR a name/alias
            
            Returns
            -------
            An integer Locus ID (LID)

        '''
        cur = self.m80.db.cursor()
        if isinstance(locus,str):
            # Handle the easy case where we have a name
            result = cur.execute(
                ' SELECT LID FROM loci WHERE name = ?',
                (locus,)
            ).fetchone()
            if result is None:
                raise MissingLocusError
            else:
                LID = result[0]
        else:
            # Try to get the LID by using the hash
            possible_lids = [x[0] for x in \
                cur.execute(
                    'SELECT LID FROM loci WHERE hash = ?',(hash(locus),)
                ).fetchall()
            ]
            # If the hash is not in the db, there is no LID
            if len(possible_lids) == 0:
                raise MissingLocusError
            # If there is one possible hash, return the LID
            elif len(possible_lids) == 1:
                LID = possible_lids[0]
            # Iterate through the loci and find the right one
            else: #pragma: no cover
                # TODO: Write a test for this case
                LID = None
                for lid in possible_lids:
                    loc = self._get_locus_by_LID(lid)
                    if locus == loc:
                        LID = lid
                        break
                if LID is None: 
                    raise MissingLocusError
        return LID

    @invalidates_primary_loci_cache
    def add_locus(self, locus: Locus, cur=None, primary_type='gene') -> int:
        '''
            Add locus to the database. 

            Parameters
            ----------
            locus : a Locus object
                This locus will be added to the db
            cur : a db cursor
                An optional cursor object to use. If none, a 
                cursor will be created. If importing many 
                loci in a loop, use a bulk transaction.

            Returns
            -------
            The locus ID (LID) of the freshly added locus
        '''
        if cur is None:
            cur = self.m80.db.cursor()
        # insert the core feature data
        core,attrs = locus.as_record()
        cur.execute(
            '''
            INSERT INTO loci 
                (chromosome,start,end,source,feature_type,strand,frame,name,hash)
                VALUES (?,?,?,?,?,?,?,?,?)
            ''',
            core,
        )
        # get the fresh LID
        (LID,) = cur.execute('SELECT last_insert_rowid()').fetchone()
        if LID is None: #pragma: no cover
            # I dont know hwen this would happen without another exception being thrown
            raise ValueError(f"{locus} was not assigned a valid LID!")
        # Add the key val pairs
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
        if locus.feature_type == primary_type:
            cur.execute('''
                INSERT INTO primary_loci (LID) VALUES (?)
            ''',(LID,))
        # get childrem LIDS using recursion
        children_LIDs = [self.add_locus(c,cur=cur) for c in locus.subloci]
        # Add relationships  
        for CID in children_LIDs:
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
            (LID,locus.start,locus.end,locus.chromosome)
        )
        return LID

    def import_gff(
        self, 
        filename: str, 
        feature_type="*", 
        ID_attr="ID", 
        parent_attr='Parent',
        attr_split="="
    ) -> None:
        '''
            Imports RefLoci from a gff (General Feature Format) file.
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
        with self.m80.db.bulk_transaction() as cur:
            for l in loci:
                self.add_locus(l,cur=cur)
        self.log.info('Done!')
        return None

    def __contains__(self, locus: Locus) -> bool:
        '''
            Returns True or False based on whether or not the Locus is
            in the database. 

            Parameters
            ----------
            locus : Locus object or str (alias)
                An input locus or the name of a locus for which 
                to look up. 

            Returns
            -------
            True or False
        '''
        try:
            # If we can get an LID, it exists
            LID = self._get_LID(locus)
            return True
        except MissingLocusError as e:
            return False

    def __getitem__(self, item):
        '''
            A convenience method to extract loci from the reference genome.
        '''
        LID = self._get_LID(item)
        return self._get_locus_by_LID(LID)

    def __iter__(self):
        return (self._get_locus_by_LID(l) for l in self._primary_LIDS())

    @contextmanager
    def filter_feature_type(self,feature_type):
        '''
        Switch primary feature types within a context block.
        Useful for temporary switching in a with block.
        
        >>> ref = locuspocus.RefLoci('test')
        >>> with ref.filter_feature_type('exon'):
                ...
                ... 'Do something with exons here'
                ...
        >>> 'primary type is now back to whatever it was before' 
        '''
        cur_lids = self._primary_LIDS()
        self.set_primary_feature_type(feature_type,clear_previous=True)
        yield
        with self.m80.db.bulk_transaction() as cur:
            cur.execute('DELETE FROM primary_loci')
            cur.executemany('''
                INSERT INTO primary_loci (LID) 
                VALUES (?)
                ''', ((x,) for x in cur_lids)
            )
        # reset the LID cache here
        self._primary_LIDS.cache_clear()


    @invalidates_primary_loci_cache
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
        with self.m80.db.bulk_transaction() as cur: 
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
        the RefLoci database
        
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
        cur = self.m80.db.cursor()
        primary_ftypes = [x[0] for x in cur.execute('''
            SELECT DISTINCT feature_type 
            FROM primary_loci p 
                JOIN loci l ON p.LID = l.LID;
        ''').fetchall()]
        ndict = dict()
        root = Node(self.name)
        for n in primary_ftypes:
            ndict[n] = Node(n,parent=root)
        ftypes = cur.execute('''
            SELECT DISTINCT p.feature_type,c.feature_type 
            FROM relationships r 
                JOIN loci p ON r.parent = p.LID 
                JOIN loci c ON r.child = c.LID;
        ''').fetchall()
        # Create Nodes
        for p,c in ftypes:
            ndict[c] = Node(c)
        for p,c in ftypes:
            if p in ndict:
                ndict[c].parent = ndict[p]
            else:
                ndict[c].parent = root
        if print_tree is True:
            print(RenderTree(root))
        return root

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
        cur = self.m80.db.cursor()
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
        cur = self.m80.db.cursor()
        LIDS = cur.execute('''
            SELECT l.LId FROM positions p, primary_loci l
            WHERE p.chromosome = ?
            AND p.start < ?
            AND p.end > ?
            AND p.LID = l.LID
        ''',(locus.chromosome,locus.start,locus.end))
        for x, in LIDS:
            yield self._get_locus_by_LID(x) 

    @invalidates_primary_loci_cache
    def _nuke_tables(self):
        cur = self.m80.db.cursor()
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

    def _initialize_tables(self):
        '''
        Initializes the Tables holding all the information
        about the Loci.
        '''
        cur = self.m80.db.cursor()
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
                name TEXT,
                hash INT
                
            );
            '''
        )

        cur.execute('''
            CREATE INDEX IF NOT EXISTS locus_id ON loci (name);
            CREATE INDEX IF NOT EXISTS locus_chromosome ON loci (chromosome);
            CREATE INDEX IF NOT EXISTS locus_start ON loci (start);
            CREATE INDEX IF NOT EXISTS locus_feature_type ON loci (feature_type);
            CREATE INDEX IF NOT EXISTS locus_end ON loci (end);
            CREATE INDEX IF NOT EXISTS locus_hash ON LOCI (hash);
            '''
        )
        cur.execute(
        # Create a table that contains loci attribute mapping
        '''
            CREATE TABLE IF NOT EXISTS loci_attrs (
                LID INT NOT NULL,
                key TEXT,
                val TEXT,
                FOREIGN KEY(LID) REFERENCES loci(LID),
                UNIQUE(LID,key)
            );
            '''
        )
        cur.execute(
        # Create a table with parent-child relationships        
        '''
            CREATE TABLE IF NOT EXISTS relationships (
                parent INT,
                child INT,
                FOREIGN KEY(parent) REFERENCES loci(LID),
                FOREIGN KEY(child) REFERENCES loic(LID)
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
                LID INT PRIMARY KEY,
                FOREIGN KEY (LID) REFERENCES loci(LID)
            );
        '''
        )



