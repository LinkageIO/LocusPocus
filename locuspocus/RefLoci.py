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

from .Locus import Locus
from .Exceptions import ZeroWindowError

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


class RefLoci(Freezable):
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
    def _primary_LIDS(self):
        '''
        A cached list of primary Locus IDs (LIDS) available in the RefLoci database.
        This list can change if the primary locus changes, for example through
        the `RefLoci.set_primary_feature_type()` method. 
        
        Note: this method takes no arguments.
        '''
        self.log.info("Caching LIDS from database")
        LIDS = [x[0] for x in self._db.cursor().execute('SELECT LID FROM primary_loci')]
        return LIDS


    def __len__(self):
        '''
        Returns the number of primary loci in the dataset. Gets called by the len()
        built-in function

        >>> len(ref)
        42
        '''
        return len(self._primary_LIDS())

    def _get_locus_by_LID(self,LID: int) -> Locus:
        '''
        Get a locus by its LID
        
        Parameters
        ----------
        LID : int
            A Locus ID. These are assigned to Locus objects when
            they are added to the RefLoci database.

        Returns
        -------
        The Locus object corresponding to the LID.

        Raises
        ------
        `TypeError` if there is no Locus in the database with that LID.
        '''
        try:
            cur = self._db.cursor()
            chromosome,start,end,source,feature_type,strand,frame,name = cur.execute('''
                SELECT chromosome, start, end, source, feature_type, strand, frame, name 
                FROM loci  
                WHERE LID = ?
            ''',(LID,)).fetchone()
            children = [x[0] for x in cur.execute('''
                SELECT child FROM relationships
                WHERE parent = ?
            ''',(LID,)).fetchall()]
            locus = Locus(
                chromosome=chromosome,
                start=start,
                end=end,
                source=source,
                feature_type=feature_type,
                strand=strand,
                frame=frame,
                subloci=children,
                name=name,
                refloci=self,
                _LID=LID
            )
            return locus
        except TypeError as e:
            raise ValueError(f'There is no locus in the database with LID of {LID}')

    def _get_LID(self,locus: Locus) -> int:
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
        cur = self._db.cursor()
        if isinstance(locus,str):
            # Handle the easy case where we have a name
            result = cur.execute(
                ' SELECT LID FROM loci WHERE name = ?',
                (locus,)
            ).fetchone()
            if result is None:
                raise ValueError('Locus not in database!')
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
                raise ValueError('Locus not in database!')
            # If there is one possible hash, return the LID
            elif len(possible_lids) == 1:
                LID = possible_lids[0]
            # Iterate through the loci and find the right one
            else:
                LID = None
                for lid in possible_lids:
                    loc = self._get_locus_by_LID(lid)
                    if locus == loc:
                        LID = lid
                        break
                if LID is None:
                    raise ValueError(f"NO LID for locus {locus}")
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
            cur = self._db.cursor()

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
        if LID is None:
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
        with self._bulk_transaction() as cur:
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
        except ValueError as e:
            return False

    def __getitem__(self, item):
        '''
            A convenience method to extract loci from the reference genome.
        '''
        LID = self._get_LID(item)
        return self._get_locus_by_LID(LID)

    def __iter__(self):
        return (self._get_locus_by_LID(l) for l in self._primary_LIDS())


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
            distinct : bool (default: True)
                If True, the return set will not contain duplicates
            primary_only : bool (default: True)
                If True, the return set will not include sub-features
                such as exons. Typically, primary features include genes.
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
        while len(loci) < n:
            try:
                rand_LID = random.choice(LIDS) 
                loci.add(self._get_locus_by_LID(rand_LID))
            except ValueError as e:
                self.log.info(f'hit a tombstone! LID {rand_LID}')
                continue
        loci = list(loci)
        if autopop and len(loci) == 1:
            loci = loci[0]
        return loci

    def feature_types(self, print_tree=True):
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
        cur = self._db.cursor()
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

    #@accepts_loci 
    def within(self, locus, partial=False):
        '''
        Returns the Loci that are within the input locus.
        By default, the entire locus needs to be within the
        coordinates of the input locus, however, this can be
        toggled with the `partial` argument.
        
        NOTE: this loci are returned in order of 3' to 5'
              based on the strand of the input locus

        partial=False (default):
              start           end
        -------[***************]-------------
            nnnnn yyyy yyy   nnnnnnnnnnn    

        partial=True:
              start           end
        -------[***************]-------------
            yyyyy yyyy yyy   yyyyyyyyyyy    

        Parameters
        ----------
        locus : Locus
            A locus defining the genomic interval to 
            extract loci within.
        partial : bool
            When True, include loci that partially overlap.
            See example above.
        '''
        # set up variables to use based on 'partial' flag
        cur = self._db.cursor()
        # Calculate the correct strand orientation
        if locus.strand == '+':
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
            raise ValueError(f'The input locus had an invalid strand: {locus.strand}')
        query = f'''
            SELECT l.LID FROM primary_loci p, loci l 
            INDEXED BY {index} 
            WHERE l.chromosome = '{locus.chromosome}' 
            AND {anchor} 
            AND l.LID = p.LID 
            ORDER BY {order}; 
        '''
        LIDS = cur.execute(query)
        # if the query is empty, yield nothing  
        if LIDS is None:
            yield
        else:
            for x, in LIDS:
                yield self._get_locus_by_LID(x) 

       
    @accepts_loci
    def upstream_loci(
        self, locus, n=np.inf, max_distance=np.inf, partial=False
    ):
        '''
            Find loci upstream of a locus.

            Loci are ordered so that the nearest loci are
            at the beginning of the list.

            NOTE: This respects the strand of the locus

            partial=False
            nnn  nnnnnnn   nnnnnn   yyyyy  nnnnnn nnnn nnnn nnnnnnnn
            partial=True
            nnn  nnnnnnn   yyyyyy   yyyyy  yyyyyy nnnn nnnn nnnnnnnn
                                            
                                           start             end
                -----------------------------x++++++++++++++++x--
                             <_______________| Window (upstream)


            Parameters
            ----------
                locus : Locus object
                    The locus object for which to fetch the upstream
                    loci
                max_distance : float (default=np.inf)
                    The maximum distance 
        '''
        # calculate the start and stop anchors 
        start,end = sorted([locus.start, locus.upstream(max_distance)])
        # create a dummy locus
        upstream_region = Locus(locus.chromosome,start,end)
        # return loci within dummy locus coordinates
        loci = self.within(upstream_region, partial=partial)
        for i,x in enumerate(loci,start=1):
            if i > n:
                break
            yield x
        #return (x for i,x in enumerate(loci) if i < n)

    @accepts_loci
    def downstream_loci(
        self, locus, n=np.inf, max_distance=np.inf, partial=False
    ):
        '''
            Returns loci downstream of a locus. 
            
            Loci are ordered so that the nearest loci are 
            at the beginning of the list.

            NOTE: this respects the strand of the locus
            
            partial=False (default)
            nnn  nnnnnnn   nnn   nnnnn  yyyy  yyyyyy yyyy nnnnnn  nnnnn
            partial=True
            nnn  nnnnnnn   nnn   yyyyy  yyyy  yyyyyy yyyy yyyyyy  nnnnn
               start             end
              ---x****************x--------------------------------
                                  |________________________^ Window (downstream)
        '''
        # calculate the start and stop anchors 
        start,end = sorted([locus.end, locus.downstream(max_distance)])
        # create a dummy locus
        upstream_region = Locus(locus.chromosome,start,end)
        # return loci within dummy locus coordinates
        loci = self.within(upstream_region, partial=partial)
        for i,x in enumerate(loci,start=1):
            if i > n:
                break
            yield x

    def flanking_loci(
        self,
        locus,
        n=np.inf,
        max_distance=np.inf,
        partial=False
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
        return (
            self.upstream_loci(locus,n=n,max_distance=max_distance,partial=partial),
            self.downstream_loci(locus,n=n,max_distance=max_distance,partial=partial)
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
        if LIDS is None:
            loci = []
        else:
            loci = [self._get_locus_by_LID(x) for (x,) in LIDS]
        return loci

    def candidate_loci(
        self,
        locus,
        nflank=2,
        max_distance=np.inf,
        chain=True,
        # Advanced options
        include_parent_locus=False,
        include_parent_attrs=False,
        include_num_intervening=False,
        include_rank_intervening=False,
        include_num_siblings=False,
        include_SNP_distance=False,
        attrs=None,
        return_table=False
    ):
        '''
            Locus to locus mapping. This method accepts an iterable of locus 
            objects (loci) and 

            Return loci between locus start and stop, plus additional
            flanking loci within window (up to flank_limit).

            Parameters
            ----------
            loci : locuspocus.Locus (also handles an iterable containing Loci)
                a locus or iterable of loci
            flank_limit : int (default : 2)
                The total number of flanking loci **on each side**
                considered a candidate surrounding a locus
            chain : bool (default : true)
                Calls itertools chain on results before returning
            window_size : int (default: None)
                Optional parameter used to extend or shorten a locus
                window from which to choose candidates from. If None,
                the function will resort to what is available in the
                window attribute of the Locus.
            include_parent_locus : bool (default: False)
                Optional parameter which will update candidate loci
                'attr' attribute with the id of the parent locus
                which contains it.
            include_parent_attrs : iterable (default: False)
                Optional parameter to include attributes from the parent
                locus. Parent locus attrs specified here will be included.
                If effective loci is > 1, the maximum value will be
                included. E.g. - including the SNP effect size with
                candidate loci.
            include_num_intervening : bool (default: False)
                Optional argument which adds an attribute to each
                candidate loci containing the rank of each loci
                as a function of distance away from the parent
                locus. (i.e. the closest candidate is 1 and the
                furthest candidate is n)
            include_rank_intervening : bool (default: False)
                Optional argument which adds the rank of each
                candidatea as a funtion of distance from the parent
                Locus. i.e. The closest loci is ranked 1 and the furthest
                loci is ranked n.
            include_num_siblings : bool (default: False)
                Optional argument which adds an attribute to each
                candidate loci containing the number of total
                candidates (siblings) identifies at the locus.
            include_SNP_distance : bool (default:False)
                Include the distance from the canadidate loci and
                the parent SNP
            attrs : dict (default: None)
                An optional dictionary which will be updated to each
                candidate loci attr value.
            return_table : bool(default: False)
                If True, return a Pandas table (DataFrame)

            Returns
            -------
            a list of candidate loci (or list of lists if chain is False)

        '''
        pass
       ## make sure to convert generators to lists
       #within = list(self.loci_within(locus,partial=False))
       #up, down = map(list,self.flanking_loci(
       #    locus, n=flank_limit, partial=True
       #))

       ## This always returns candidates together, if
       ## you want specific up,within and down loci
       ## use the specific methods
       #candidates = sorted(up_loci+loci_within+down_loci))
       #return candidates
       ## include the number of effective loci
       #if include_rank_intervening == True:
       #    ranks = sp.stats.rankdata(
       #        [abs(x.center_distance(locus)) for x in candidates]
       #    )
       ## Iterate through candidate loci and propagate the
       ## parental info
       #for i, cand in enumerate(candidates):
       #    # include parent locus id if thats specified
       #    if include_parent_locus == True:
       #        cand.update({"parent_locus": locus.id})
       #    if include_rank_intervening == True:
       #        cand.update({"intervening_rank": ranks[i]})
       #    # update all the parent_attrs
       #    if include_parent_attrs and len(include_parent_attrs) > 0:
       #        if "all" in include_parent_attrs:
       #            include_parent_attrs = locus.attr.keys()
       #        for attr in include_parent_attrs:
       #            attr_name = "parent_{}".format(attr)
       #            cand.update({attr_name: locus[attr]})
       #if include_num_intervening == True:
       #    num_down = 0
       #    num_up = 0
       #    # Sort the loci by their distance from the locus
       #    cands_with_distances = [
       #        (cand, abs(cand.center_distance(locus))) for cand in candidates
       #    ]
       #    cands_with_distances = sorted(cands_with_distances, key=lambda x: x[1])
       #    for cand, distance in cands_with_distances:
       #        if locus.within(cand):
       #            cand.update({"num_intervening": -1})
       #        elif cand.center >= locus.center:
       #            cand.update({"num_intervening": num_down})
       #            num_down += 1
       #        elif cand.center <= locus.center:
       #            cand.update({"num_intervening": num_up})
       #            num_up += 1
       #if include_num_siblings == True:
       #    for cand in candidates:
       #        cand.update({"num_siblings": len(candidates)})
       #if include_SNP_distance == True:
       #    for cand in candidates:
       #        distance = abs(cand.center_distance(locus))
       #        cand.update({"SNP_distance": distance})
       #if attrs is not None:
       #    for cand in candidates:
       #        cand.update(attrs)
       #if return_table == True:
       #    candidates = pd.DataFrame([x.as_dict() for x in candidates])
       #return candidates

    def _rtree_within(self, locus, partial=False):
        '''
            Implements the wihtin function using the R*Tree index
            This method has the same functionality as RefLoci.within()
            but just used a different SQLite approach. The within()
            method appears to be faster, but I want to keep this
            around in the commit history for a bit
        '''
        cur = self._db.cursor()
        if partial == False:
            LIDS = cur.execute('''
                SELECT l.LID FROM positions p, primary_loci l 
                WHERE p.chromosome = ? 
                AND p.start > ? 
                AND p.end < ?
                AND p.LID = l.LID
                ORDER BY p.start;
            ''',(locus.chromosome,locus.start,locus.end))
        else:
            LIDS = cur.execute('''
                SELECT l.LID FROM positions p, primary_loci l 
                WHERE p.chromosome = ?
                AND p.start < ? 
                AND p.end > ?
                AND p.LID = l.LID
                ORDER BY p.start;
            ''',(locus.chromosome,locus.end,locus.start))
        yield from (self._get_locus_by_LID(x) for (x,) in LIDS)
 

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

    def _initialize_tables(self):
        '''
        Initializes the Tables holding all the information
        about the Loci.
        '''
        cur = self._db.cursor()
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



