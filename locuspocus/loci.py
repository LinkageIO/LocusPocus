#!/usr/bin/python3
import apsw
import gzip
import random
import logging

import numpy as np
import minus80 as m80

from typing import (
    List, Optional
)

from pathlib import Path
from minus80 import Freezable
from functools import wraps

from .locus import Locus
from .locusview import LocusView
from .exceptions import MissingLocusError,StrandError

__all__ = ['Loci']

log = logging.getLogger(__name__)

# --------------------------------------------------
#       Decorators
# --------------------------------------------------

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


# --------------------------------------------------
#       Class Definition
# --------------------------------------------------

class Loci(Freezable):
    '''
        Loci are more than the sum of their parts. They have a name and
        represent something bigger than theirselves. They are important. They
        live on the disk in a database.
    '''

    def __init__(
            self, 
            name: str, 
            rootdir: Optional[str] = None
        ):
        '''
            Initialize a new Locus object

            Parameters
            ----------

            name : str 
                The name of the now Loci object
            rootdir : str
                The base directory to store the files related to the dataset
                If not specified, the default will be taken from the config file
        '''
        # set up the freezable API
        super().__init__(name, rootdir=rootdir)
        self.name = name
        self._initialize_tables()
        self._cached_LIDs = None

    @property
    def _LIDs(self) -> List[int]:
        if self._cached_LIDs is None:
            self._cached_LIDs = [LID for (LID,) in self.m80.db.cursor().execute(
                'SELECT LID FROM loci'    
            )]
        return self._cached_LIDs

    def __len__(self) -> int:
        '''
        Returns the number of loci in the dataset. 
        Gets called by the len().

        >>> len(ref)
        42
        '''
        (num_rows,) = self.m80.db.cursor().execute(
            'SELECT COUNT(*) FROM LOCI;'        
        ).fetchone()
        return num_rows

    def _get_locus_by_LID(self,LID: int) -> LocusView:
        '''
        Get a locus by its LID
        
        Parameters
        ----------
        LID : int
            A Locus ID. These are assigned to Locus objects when
            they are added to the Loci database.

        Returns
        -------
        The LocusView.

        Raises
        ------
        `MissingLocusError` if there is no Locus in the database with that LID.
        '''
        lid_exists, = self.m80.db.cursor().execute(
            'SELECT COUNT(*) FROM loci WHERE LID = ? ',
            (LID,)
        ).fetchone()
        if lid_exists == 0:
            raise MissingLocusError(f"Cannot find Locus for LID: {LID}")
        return LocusView(LID,self)

    def _get_LID(self,locus: Locus) -> int: #pragma: no cover
        '''
            Return the Locus Identifier used internally by Loci

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
                raise MissingLocusError(f"Cannot find LID for Locus: {locus}")
            else:
                LID = result[0]
        else:
            raise MissingLocusError(f"Cannot find LID for Locus: {locus}")
        return LID

    def add_locus(
        self, 
        locus: Locus, 
        cur=None, 
    ) -> int:
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
                (chromosome,start,end,source,feature_type,strand,frame,name)
                VALUES (?,?,?,?,?,?,?,?)
            ''',
            core,
        )
        # get the locus LID
        (LID,) = cur.execute('SELECT last_insert_rowid()').fetchone()

        if LID is None: # pragma: no cover
            # I dont know when this would happen without another exception being thrown
            raise ValueError(f"{locus} was not assigned a valid LID!")
        # Add the attrs
        for key,val in attrs.items():
            cur.execute('''
                INSERT INTO loci_attrs 
                    (LID,key,val) 
                    VALUES (?,?,?)
                ''', 
                (LID,key,val)
            )
        # Add subloci information
        self._add_subloci(
            root_LID=LID,
            parent_LID=None,
            subloci=locus.subloci,
            cur=cur
        )
        # Add the position to the R*Tree
        cur.execute(
            '''
            INSERT INTO positions (LID,start,end,chromosome) VALUES (?,?,?,?)
            ''',
            (LID,locus.start,locus.end,locus.chromosome)
        )
        return LID

    def _add_subloci(
        self,
        root_LID: int,
        parent_LID: int,
        subloci: "Loci",
        cur: Optional[apsw.Connection.cursor] = None
    ):
        if subloci.empty:
            return None
        for l in subloci:
            # Extract info
            core,attrs = l.as_record()
            # add the sublocus
            cur.execute('''
                INSERT INTO subloci
                    (root_LID,parent_LID,chromosome,start,end,source,feature_type,strand,frame,name)
                    VALUES (?,?,?,?,?,?,?,?,?,?)''',
                (root_LID,parent_LID)+core 
            )
            (LID,) = cur.execute('SELECT last_insert_rowid()').fetchone()
            # add the attrs
            for key,val in l.attrs.items():
                cur.execute('''
                    INSERT INTO subloci_attrs 
                    (LID,key,val) 
                    VALUES (?,?,?)
                    ''', 
                    (LID,key,val)
                )
            # recurse with updated parentLID
            self._add_subloci(
                root_LID=root_LID,
                parent_LID=LID,
                subloci=l.subloci,
                cur=cur
            )

    def import_gff(
        self, 
        filename: str, 
        /,
        ID_attr: str = "ID", 
        parent_attr: str = 'Parent',
        attr_split: str = "=",
        skip_feature_types: Optional[List[str]] = None
    ) -> None:
        '''
            Imports Loci from a gff (General Feature Format) file
            to the current Loci object.

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
            skip_feature_types : Optional[List[str]]
                Optionally, provide a list of feature_types to skip during the
                import. For instance, some GFFs will provide features for 
                Chromosomes, which can lead to strange behaviors.
        '''
        log.info(f"Importing Loci from {filename}")
        if filename.endswith(".gz"):
            IN = gzip.open(filename, "rt")
        else:
            IN = open(filename, "r")
        loci = []
        current_locus = None
        total_loci = 0
        for i,line in enumerate(IN):
            total_loci += 1
            # skip comment lines
            if line.startswith("#"):
                continue
            locus = Locus.from_gff_line(line, ID_attr=ID_attr, parent_attr=parent_attr, attr_split=attr_split)
            # Check to see if we are in a top level locus
            if skip_feature_types and locus.feature_type in skip_feature_types:
                continue
            if parent_attr not in locus.attrs:
                current_locus = locus
                loci.append(locus)
            else:
                # add the sublocus to the current locus
                current_locus.add_sublocus(locus,find_parent=True)
        log.info((
            f'Found {len(loci)} loci, adding to database'
        ))
        IN.close()
        with self.m80.db.bulk_transaction() as cur:
            for l in loci:
                self.add_locus(l,cur=cur)
        log.info('Done!')
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
            self._get_LID(locus)
            return True
        except MissingLocusError:
            return False

    def __getitem__(self, item):
        '''
            A convenience method to extract a locus from 
            the loci.
        '''
        LID = self._get_LID(item)
        return self._get_locus_by_LID(LID)

    def __iter__(self):
        LIDs = self.m80.db.cursor().execute('''
            SELECT LID from loci
        ''')
        return (self._get_locus_by_LID(l) for (l,) in LIDs)

    def rand(
        self, 
        n: int = 1, 
        distinct: bool = True, 
        autopop: bool = True
    ):
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
        loci = set()
        if n > len(self._LIDs):
            raise ValueError('More than the maximum loci in the database was requested')
        if distinct == True:
            LIDs = random.sample(self._LIDs,n)
        else:
            LIDs = random.choices(self._LIDs,n)
        loci = [self._get_locus_by_LID(x) for x in LIDs]
        if autopop and len(loci) == 1:
            loci = loci[0]
        return loci
    
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
            SELECT l.LID FROM loci l 
            INDEXED BY {index} 
            WHERE l.chromosome = '{locus.chromosome}' 
            AND {anchor} 
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
            SELECT LID FROM positions
            WHERE chromosome = ?
            AND start < ?
            AND end > ?
        ''',(locus.chromosome,locus.start,locus.end))
        for x, in LIDS:
            yield self._get_locus_by_LID(x) 

    def _nuke_tables(self):
        cur = self.m80.db.cursor()
        cur.execute(
            '''
                DROP TABLE IF EXISTS loci;
                DROP TABLE IF EXISTS subloci;
                DROP TABLE IF EXISTS loci_attrs;
                DROP TABLE IF EXISTS subloci_attrs;
                DROP TABLE IF EXISTS positions;
            '''
        )
        self._initialize_tables()

    def _initialize_tables(self):
        '''
            Initializes the Tables holding all the information
            about the Loci.
        '''
        cur = self.m80.db.cursor()
        cur.execute('''
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

                name TEXT
            );
            CREATE INDEX IF NOT EXISTS locus_LID on loci (LID);
            CREATE INDEX IF NOT EXISTS locus_id ON loci (name);
            CREATE INDEX IF NOT EXISTS locus_chromosome ON loci (chromosome);
            CREATE INDEX IF NOT EXISTS locus_start ON loci (start);
            CREATE INDEX IF NOT EXISTS locus_end ON loci (end);
            CREATE INDEX IF NOT EXISTS locus_feature_type ON loci (feature_type);
        ''')

        cur.execute('''
             CREATE TABLE IF NOT EXISTS subloci (
                LID INTEGER PRIMARY KEY AUTOINCREMENT,
                
                /* Store the locus values  */
                root_LID INTEGER,
                parent_LID INTEGER,

                chromosome TEXT NOT NULL,
                start INTEGER NOT NULL,
                end INTEGER,

                source TEXT,
                feature_type TEXT,
                strand TEXT,
                frame INT,

                name TEXT 
            );
            CREATE INDEX IF NOT EXISTS subloci_LID ON subloci (LID);
            CREATE INDEX IF NOT EXISTS subloci_root_LID ON subloci (root_LID);
            CREATE INDEX IF NOT EXISTS subloci_parent_LID ON subloci (parent_LID);
        ''')

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
            CREATE INDEX IF NOT EXISTS loci_attrs_LID ON loci_attrs (LID);
            CREATE INDEX IF NOT EXISTS loci_attrs_LID_key ON loci_attrs (LID,key);
            '''
        )

        cur.execute(
        # Create a table that contains loci attribute mapping
        '''
            CREATE TABLE IF NOT EXISTS subloci_attrs (
                LID INT NOT NULL,
                key TEXT,
                val TEXT,
                FOREIGN KEY(LID) REFERENCES subloci(LID),
                UNIQUE(LID,key)
            );
            CREATE INDEX IF NOT EXISTS subloci_attrs_LID ON subloci_attrs (LID);
            CREATE INDEX IF NOT EXISTS subloci_attrs_LID_key ON subloci_attrs (LID,key);
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

    # --------------------------------------------------
    #       factory methods
    # --------------------------------------------------
    @classmethod
    def from_gff(
        cls,
        name: str,
        gff_file: str, 
        /,
        rootdir: Optional[str] = None,
        ID_attr: str = "ID", 
        parent_attr: str = 'Parent',
        attr_split: str = "=",
        overwrite: bool = False,
        skip_feature_types: Optional[List[str]] = None
    ) -> "Loci":
        '''
            Create a new Loci object from a GFF file.

            Parameters
            ----------
            name : str
                The name of the resultant Loci object
            gff_file: str
                The path to the gff file
            ID_attr : str (default: ID)
                The key in the attribute column which designates the ID or
                name of the feature.
            parent_attr : str (default: Parent)
                The key in the attribute column which designates the Parent of
                the Locus
            attr_split : str (default: '=')
                The delimiter for keys and values in the attribute column
            overwrite : bool (default: False)
                If True, any existing dataset with this name will be
                deleted prior to building.
            skip_features_types : Optional[List[str]]
                Optionally, provide a list of feature_types to skip during the
                import. For instance, some GFFs will provide features for 
                Chromosomes, which can lead to strange behaviors.
        '''
        gff_file = Path(gff_file)
        if overwrite:
            m80.Tools.delete('Loci',name)
        # Do some checks
        if m80.exists('Loci',name):
            raise ValueError(
                f'Loci.{name} exists. Cannot use factory '
                f'methods on existing datasets.'
            )
        if not gff_file.exists():
            raise FileNotFoundError(f"{gff_file} does not exist")
        # Create and import the features
        loci = cls(name, rootdir=rootdir)
        loci.import_gff(
            str(gff_file),
            ID_attr=ID_attr,
            parent_attr=parent_attr,
            attr_split=attr_split,
            skip_feature_types=skip_feature_types
        )
        return loci

    @classmethod
    def from_loci(
        cls,
        name: str,
        source_loci: "Loci",
        /,
        rootdir: Optional[str] = None,

    ) -> "Loci":
        '''
            Efficiently create a new Loci object based on a 
            list of filtered Locus names and a source Loci object. 

            Parameters
            ----------
            name : str 
                The name of the now Loci object
            source_loci : locuspocus.Loci
                The Loci object to filter from
            names : List[str]
                A list of names to filter from the source Loci
            rootdir : Optional[str]
                
        '''
        # Do some checks
        if m80.exists('Loci',name):
            raise ValueError(
                f'Loci.{name} exists. Cannot use factory '
                f'methods on existing datasets.'
            )

        loci = cls(name, rootdir=rootdir)
   
        with loci.m80.db.bulk_transaction() as cur:
            for l in source_loci:
                loci.add_locus(l,cur=cur)
        return loci

