#!/usr/bin/python3
import itertools
import collections
import random
import numpy as np
import scipy as sp
import math
import gzip
import re
import logging 

from minus80 import Freezable
from collections import defaultdict
from functools import lru_cache,wraps

from .LocusDist import LocusDist
from .Locus import Locus
from .Exceptions import ZeroWindowError

def memoize(obj):
    cache = obj.cache = {}
    @wraps(obj)
    def memoizer(*args, **kwargs):
        # Give us a way to clear the cache
        if 'clear_cache' in kwargs:
            cache.clear()
        # This wraps the calling of the memoized object
        key = str(args) + str(kwargs)
        if key not in cache:
            cache[key] = obj(*args, **kwargs)
        return cache[key]
    return memoizer


class Loci(Freezable):
    '''
        Just a bunch of Locuses. Loci are more than the sum of their 
        parts. They have a name and represent something bigger than
        theirselves. They are important. They live on the disk in a
        database.
    '''
    log = logging.getLogger(__name__)
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
                    '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
                )
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(logging.INFO)

    def __init__(self,name):
        super().__init__(name)
        self.name = name
        self._initialize_tables()
    def __len__(self):
        '''
            Returns the number of Loci in the database
        '''
        return self.num_loci()
    def __contains__(self,obj):
        ''' flexible on what you pass into the 'in' function '''
        if isinstance(obj,Locus):
            # you can pass in a Locus object (this expression
            # should ALWAYS be true if you
            # created Locus object from this Loci object)
            if self._db.cursor().execute(
                '''SELECT COUNT(*) FROM loci WHERE id = ?''',
                (obj.id,)).fetchone()[0] == 1:
                return True
            else:
                return False
        elif isinstance(obj,str):
            # Can be a string object
            if self._db.cursor().execute('''
                SELECT COUNT(*) FROM loci WHERE id = ?''',
                (obj,)).fetchone()[0] == 1:
                return True
            elif self._db.cursor().execute(''' 
                SELECT COUNT(*) FROM aliases WHERE alias = ?''',
                (obj,)).fetchone()[0] == 1:
                return True
            else:
                return False
        else:
            raise TypeError('Cannot test for containment for {}'.format(obj))
    def __getitem__(self,item):
        '''
            A convenience method to extract loci from the reference genome.
        '''
        # NOTE: Dont LRU cache this, it gets cached in from_id
        if isinstance(item,str):
            return self.get_locus_from_id(item)
        # Allow for iterables of items to be passed in
        else:
            return self.get_loci_from_ids(item)

    def add_loci(self,loci):
        '''
            Add loci to the database

            Parameters
            ----------
            loci : an iterable of loci
                These loci are added to the database

            Returns
            -------
            None
        '''
        if isinstance(loci,Locus):
            # Its actually just one locus
            locus = loci
            self._db.cursor().execute('''
            INSERT OR REPLACE INTO loci VALUES (?,?,?,?)
            ''',(locus.name, locus.chrom, locus.start, locus.end))
            self._db.cursor().executemany('''
            INSERT OR REPLACE INTO loci_attrs VALUES (?,?,?)
            ''',[(locus.id,key,val) for key,val in locus.attr.items()])
        else:
            try:
                # support adding lists of loci
                self.log.info('Adding {} Genes info to database'.format(len(loci)))
                cur = self._db.cursor()
                cur.execute('BEGIN TRANSACTION')
                cur.executemany(
                    'INSERT OR REPLACE INTO loci VALUES (?,?,?,?)',
                    ((x.name,x.chrom,x.start,x.end) for x in loci)
                )
                self.log.info('Adding Gene attr info to database')
                cur.executemany(
                    'INSERT OR REPLACE INTO loci_attrs VALUES (?,?,?)',
                    ((x.id,key,val) for x in loci for key,val in x.attr.items())
                )
                cur.execute('END TRANSACTION')
            except Exception as e:
                cur.exectue('ROLLBACK')

    def remove_locus(self,item):
        '''
            Remove loci from the database

            Parameters
            ----------
            loci : an iterable of IDs or a single ID

            Returns
            -------
            The number of loci removed
        '''
        if isinstance(item, str):
            # We have a single locus
            cur = self._db.cursor()
            cur.execute("DELETE FROM loci WHERE id = ?",(item,))
        elif isinstance(item, tuple):
            # We have coordinates
            pass
        else:
            raise ValueError('Cannot find: {}. Must be an ID or coordinate tuple.')


    def add_gff(self,filename,
                 locus_feature='gene',
                 ID_attr='ID',attr_split='='):
        '''
            Imports Loci from a gff (General Feature Format) file.
            See more about the format here:
            http://www.ensembl.org/info/website/upload/gff.html

            Parameters
            ----------

            filename : str
                The path to the GFF file.
            name : str
                The name if the RefGen object to be stored in the core
                camoco database.
            locus_feature : str (default: 'gene')
                The name of the feature (in column 2) that designates a 
                locus. These features will be the main object that the RefGen
                encompasses. 
            ID_attr : str (default: ID)
                The key in the attribute column which designates the ID or 
                name of the feature.
            attr_split : str (default: '=')
                The delimiter for keys and values in the attribute column
        '''
        loci = list()
        chroms = dict()
        if filename.endswith('.gz'):
            IN = gzip.open(filename,'rt')
        else:
            IN = open(filename,'r')
        for line in IN:
            #skip comment lines
            if line.startswith('#'):
                continue
            (chrom,source,feature,start,
             end,score,strand,frame,attributes) = line.strip().split('\t')
            attributes = dict([(field.strip().split(attr_split)) \
                    for field in attributes.strip(';').split(';')])
            if feature == locus_feature:
                loci.append(
                    Locus(
                        chrom,int(start),int(end),
                        attributes[ID_attr].upper().strip('"'),strand=strand,
                        **attributes
                    ).update(attributes)
                )
                # Keep track of seen chromosomes
                if chrom not in chroms:
                    chroms[chrom] = end
                else:
                    if end > chroms[chrom]:
                        chroms[chrom] = end
        IN.close()
        self.add_loci(loci)

    #----------- Internal Methods  -----------#


    @memoize
    def num_loci(self):
        '''
            Returns the number of loci in the dataset
        '''
        return self._db.cursor().execute(
           ''' SELECT COUNT(*) FROM loci'''
        ).fetchone()[0]

    def iter_loci(self):
        '''
            Iterates over loci in database.

            Returns
            -------
            A generator containing loci
        '''
        for x in self._db.cursor().execute('''
                SELECT chromosome,start,end,id FROM loci
            '''):
            yield Locus(*x)
    def random_locus(self,**kwargs):
        '''
            Returns a random locus within loci.
            Also allows passing of keyword arguments to Locus
            constructor method allowing for flexible generation.
            See Locus.__init__ for more details.

            Parameters
            ----------
            **kwargs : key,value pairs
                Extra parameters passed onto the locus init method.

            Returns
            -------
            A Locus object (camoco.Locus based)

        '''
        chrom,start,end,id = self._db.cursor().execute('''
            SELECT chromosome,start,end,id from loci ORDER BY RANDOM() LIMIT 1;
        ''').fetchone()
        return Locus(chrom,start,end,id=id,**kwargs
        )
    def random_loci(self,n,**kwargs):
        '''
            Return random loci from the database, without replacement. 

            Parameters
            ----------
            n : int

            **kwargs : key,value pairs
                Extra parameters passed onto the locus init method

            Returns
            -------
            An iterable containing n (unique) random loci

        '''
        loci_info = self._db.cursor().execute(
                "SELECT chromosome,start,end,id from loci ORDER BY RANDOM() LIMIT ?",
                (n,)
        )
        return set([Locus(chr,start,end=end,id=id,**kwargs) for \
            (chr,start,end,id) in loci_info])

    def intersection(self,loci):
        '''
            Return the subset of loci that are in the dataset.
            
            Parameters
            ----------
            loci : list-like of Locus

            Returns
            -------
            a list like object containing loci that
            are in the dataset.
        '''
        return [x for x in loci if x in self]


    # --------- Getters -----------#
    @lru_cache(maxsize=131072)
    def get_locus_from_id(self, locus_id):
        '''
            Returns a locus object from a string

            Parameters
            ----------
            locus_id : str
                ID of the locus you want to pull out

            Returns
            -------
            A single Locus Object

        '''
        cur = self._db.cursor()
        # The most common use case is a direct reference to an id
        info = cur.execute(
                'SELECT chromosome,start,end,id FROM loci WHERE id = ?', [locus_id]
               ).fetchone()
        if info != None:
            return Locus(*info)
        else:
            # Try to fetch an alias
            (locus_id,) = self._db.cursor().execute('''
                SELECT id FROM aliases 
                WHERE alias = ?
            ''',
                (locus_id,)
            ).fetchone()
            return self[locus_id]
            
    def get_loci_from_ids(self, locus_ids, check_shape=False):
        '''
            Returns a list of locus object from an iterable of id strings
            OR from a single locus id string.

            Parameters
            ----------
            locus_ids : str OR iterable of str
                ID(s) of the loci you want to pull out
            check_shape : bool (default: False)
                Check if you get back the same number of ids you
                pass in. If false (default), just give back what
                you find, ignoring erronous ids.

            Returns
            -------
            A list of locus objects if you pass in an iterable,
            otherwise a single locus

        '''
        if isinstance(locus_ids,str):
            import warnings
            warnings.warn(
                'Passing single values into RefGen.from_ids is deprecated. Use RefGen.from_id() '
                'or slicing syntax instead.'    
            )
            return self.get_locus_from_id(locus_ids)
        loci = []
        for id in locus_ids:
            try:
                loci.append(self.get_locus_from_id(id))
            except ValueError as e:
                if check_shape == False:
                    continue
                else:
                    raise e
        return loci


    def upstream_loci(self,locus,locus_limit=1000,window_size=None):
        '''
            Find loci that START upstream of a locus. 
            Loci are ordered so that the nearest loci are 
            at the beginning of the list.

            Return Genes that overlap with the upstream window,
            This includes partially overlapping loci, but NOT
            loci that are returned by the loci_within method. 

            Looks like: (y=yes,returned; n=no,not returned)

            nnn  yyyyyyy   yyyyyy   yyyyy  yyyyyy nnnn nnnn nnnnnnnn
                                             nnnn
                                           start             end
                -----------------------------x****************x--
                   ^_________________________| Window (upstream)
        '''
        if locus.window == 0 and window_size is None:
            raise ValueError(
                'Asking for upstream loci for {} which has a window size of 0'.format(
                    locus
                )
            )
        if window_size is not None:
            upstream = locus.start - window_size
        else:
            upstream = locus.upstream
        return [
            Locus(*x) \
            for x in self._db.cursor().execute('''
                SELECT chromosome,start,end,id FROM loci 
                INDEXED BY loci_start_end
                WHERE chromosome = ?
                AND start >= ?  -- Gene must end AFTER locus window (upstream) 
                AND start < ? -- Gene must start BEFORE locus
                ORDER BY start DESC
                LIMIT ?
            ''',(locus.chrom, upstream,locus.start, locus_limit)
        )]
    def downstream_loci(self,locus,locus_limit=1000,window_size=None):
        '''
            Returns loci downstream of a locus. Genes are ordered 
            so that the nearest loci are at the beginning of the list.

            Return Genes that overlap with the downstream window,
            This includes partially overlapping loci, but NOT
            loci that are returned by the loci_within method. 

            Looks like: (y=yes,returned; n=no,not returned)

            nnn  nnnnnnn   nnnnnn nnnn  yyyy  yyyyyy yyyy yyyyyy  nnnnn
               start             end
              ---x****************x--------------------------------
                                  |_______________________^ Window (downstream)
        '''
        if locus.window == 0 and window_size is None:
            raise ValueError(
                'Asking for upstream loci for {} with no window size',
                locus.id
            )
        if window_size is not None:
            downstream = locus.end + window_size
        else:
            downstream = locus.downstream

        return [
            Locus(*x) \
            for x in self._db.cursor().execute('''
                SELECT chromosome,start,end,id FROM loci
                INDEXED BY loci_start_end
                WHERE chromosome = ?
                AND start > ?
                AND start <= ?
                ORDER BY start ASC
                LIMIT ?
            ''',(locus.chrom, locus.end, downstream, locus_limit)
        )]
    def flanking_loci(self, loci, flank_limit=2,chain=True,window_size=None):
        '''
            Returns loci upstream and downstream from a locus
            ** done NOT include loci within locus **
        '''
        if isinstance(loci,Locus):
            # If we cant iterate, we have a single locus
            locus = loci
            if locus.window == 0 and window_size is None:
                raise ValueError(
                    'Asking for upstream loci for {} and no window size.',
                    locus.id
                )
            upstream_locus_limit = int(flank_limit)
            downstream_locus_limit = int(flank_limit)
            up_loci = self.upstream_loci(
                locus, locus_limit=upstream_locus_limit, window_size=window_size
            )
            down_loci = self.downstream_loci(
                locus, locus_limit=downstream_locus_limit, window_size=window_size
            )
            if chain:
                return list(itertools.chain(up_loci,down_loci))
            return (up_loci,down_loci)
        else:
            iterator = iter(loci)
            loci = [
                self.flanking_loci(locus,flank_limit=flank_limit,window_size=window_size)\
                for locus in iterator
            ]
            if chain:
                loci = list(itertools.chain(*loci))
            return loci

    def encompassing_loci(self,loci,chain=True):
        '''
            Returns the Loci encompassing the locus. In other words
            if a locus (e.g. a SNP) is inside of another locus, i.e. the
            start of the locus is upstream and the end of the locus
            is downstream of the locus boundaries, this method will
            return it. Not to be confused with candidate locus, which
            will return locus upstream and downstream surrounding a locus.

            Parameters
            ----------
            loci : an iterable of loci
                The method will return encompassing loci for each
                locus in ther iterable. If a single locus is passed,
                a single result will be returned.
            chain : bool (defualt: True)
                Specifies whether or not to chain results. If true
                a single list will be returned, if False, a result for
                each locus passed in will be returned.
        '''
        if isinstance(loci,Locus):
            return [
                Locus(*x) \
                for x in self._db.cursor().execute('''
                    SELECT chromosome,start,end,id FROM loci
                    WHERE chromosome = ?
                    AND start <= ? AND end >= ?
                    ''',
                    (loci.chrom,loci.start,loci.end))
            ]
        else:
            iterator = iter(loci)
            loci = [self.encompassing_loci(locus,chain=chain) for locus in iterator]
            if chain:
                loci = list(itertools.chain(*loci))
            return loci
    def loci_within(self,loci,chain=True):
        '''
            Returns the loci that START within a locus 
            start/end boundry.

            Looks like: (y=yes,returned; n=no,not returned)

            nnn  nnnnnnn  yyyyyy   yyyyy  yyyyyy yyyyyyy
                    start                        end
                -----x****************************x-------------

        '''
        if isinstance(loci,Locus):
            return [
                Locus(*x) \
                for x in self._db.cursor().execute('''
                    SELECT chromosome,start,end,id FROM loci
                    WHERE chromosome = ?
                    AND start >= ? AND start <= ?
                    ''',
                    (loci.chrom,loci.start,loci.end))
            ]
        else:
            iterator = iter(loci)
            loci = [self.loci_within(locus,chain=chain) for locus in iterator]
            if chain:
                loci = list(itertools.chain(*loci))
            return loci

    def candidate_loci(self, loci, flank_limit=2,
        chain=True, window_size=None, include_parent_locus=False,
        include_parent_attrs=False, include_num_intervening=False, 
        include_rank_intervening=False, include_num_siblings=False,
        include_SNP_distance=False,attrs=None,return_table=False):
        '''
            Locus to locus mapping.
            Return loci between locus start and stop, plus additional
            flanking loci (up to flank_limit)

            Parameters
            ----------
            loci : camoco.Locus (also handles an iterable containing Loci)
                a camoco locus or iterable of loci
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
        if isinstance(loci,Locus):
            # If not an iterator, its a single locus
            locus = loci
            loci_within = self.loci_within(locus)
            up_loci,down_loci = self.flanking_loci(
                locus, flank_limit=flank_limit, chain=False,
                window_size=window_size
            )

            # This always returns candidates together, if 
            # you want specific up,within and down loci
            # use the specific methods
            candidates = sorted(itertools.chain(up_loci,loci_within,down_loci))
            #include the number of effective loci
            if include_rank_intervening == True:
                ranks = sp.stats.rankdata([abs(x.center_distance(locus)) for x in candidates])
            # Iterate through candidate loci and propagate the 
            # parental info
            for i,cand in enumerate(candidates):
                # include parent locus id if thats specified
                if include_parent_locus == True:
                    cand.update({'parent_locus':locus.id})
                if include_rank_intervening == True:
                    cand.update({'intervening_rank':ranks[i]})
                # update all the parent_attrs
                if include_parent_attrs and len(include_parent_attrs) > 0:
                    if 'all' in include_parent_attrs: 
                        include_parent_attrs = locus.attr.keys()
                    for attr in include_parent_attrs:
                        attr_name = 'parent_{}'.format(attr)
                        cand.update({attr_name: locus[attr]})
            if include_num_intervening == True:
                num_down = 0
                num_up = 0
                # Sort the loci by their distance from the locus
                cands_with_distances = [(cand,abs(cand.center_distance(locus))) for cand in candidates]
                cands_with_distances = sorted(cands_with_distances,key=lambda x: x[1])
                for cand,distance in cands_with_distances:
                    if locus.within(cand):
                        cand.update({'num_intervening':-1})
                    elif cand.center >= locus.center:
                        cand.update({'num_intervening':num_down})
                        num_down += 1
                    elif cand.center <= locus.center:
                        cand.update({'num_intervening':num_up})
                        num_up += 1
            if include_num_siblings == True:
                for cand in candidates:
                    cand.update({'num_siblings':len(candidates)})
            if include_SNP_distance == True:
                for cand in candidates:
                    distance = abs(cand.center_distance(locus))
                    cand.update({'SNP_distance':distance})
            if attrs is not None:
                for cand in candidates:
                    cand.update(attrs)
            if return_table == True:
                candidates = pd.DataFrame([x.as_dict() for x in candidates])
            return candidates

        else:
            iterator = iter(sorted(loci))
            candidates = [
                # This is becoming a pain in the ass
                self.candidate_loci(
                    locus, 
                    flank_limit=flank_limit,
                    chain=chain,
                    window_size=window_size,
                    include_parent_locus=include_parent_locus,
                    include_parent_attrs=include_parent_attrs,
                    include_num_intervening=include_num_intervening,
                    include_rank_intervening=include_rank_intervening,
                    include_num_siblings=include_num_siblings,
                    include_SNP_distance=include_SNP_distance,
                    return_table=return_table,
                    attrs=attrs
                ) for locus in iterator
            ]
            if chain:
                if return_table:
                    candidates = pd.concat(candidates)
                else:
                    candidates = list(set(itertools.chain(*candidates)))
            return candidates
    def bootstrap_candidate_loci(self, loci, flank_limit=2,
        chain=True, window_size=None, include_parent_locus=False):
        '''
            Returns candidate loci which are random, but conserves
            total number of overall loci.

            Parameters
            ----------
            loci : camoco.Locus (also handles an iterable containing Loci)
                a camoco locus or iterable of loci
            flank_limit : int (default : 2)
                The total number of flanking loci **on each side**
                considered a candidate surrounding a locus
            chain : bool (default : true)
                Calls itertools chain on results before returning,
            include_parent_locus : bool (default: False)
                Optional parameter which will update candidate loci
                'attr' attribute with the id of the parent locus
                which contains it.

            Returns
            -------
            a list of candidate loci (or list of lists if chain is False)

        '''
        if isinstance(loci,Locus):
            # We now have a single locus
            locus = loci
            # grab the actual candidate loci
            num_candidates = len(
                self.candidate_loci(
                    locus, flank_limit=flank_limit,
                    chain=True, window_size=window_size
                )
            )
            if num_candidates == 0:
                return []
            # Snps a random loci from the genome
            random_locus = self.random_locus()
            # Snag the same number of candidates
            random_candidates = self.upstream_loci(
                random_locus, 
                locus_limit=num_candidates,
                window_size=10e100
            )
            if len(random_candidates) != num_candidates:
                # somehow we hit the end of a chromosome
                # or something, just recurse
                random_candidates = self.bootstrap_candidate_loci(
                    locus,flank_limit=flank_limit,chain=True
                )
            if include_parent_locus == True:
                for cand in random_candidates:
                    cand.update({'parent_locus':random_locus.id})
            return random_candidates
        else:
            # Sort the loci so we can collapse down
            locus_list = sorted(loci)
            seen = set()
            bootstraps = list()
            target = self.candidate_loci(
                locus_list,flank_limit=flank_limit,
                chain=False,window_size=window_size
            )
            target_accumulator = 0
            candidate_accumulator = 0
            #self.log.info('target: {}, loci: {}',len(target),len(locus_list))
            for i,(locus,targ) in enumerate(zip(locus_list,target)):
                # compare downstream of last locus to current locus
                candidates = self.bootstrap_candidate_loci(
                    locus, flank_limit=flank_limit, 
                    chain=True, window_size=window_size,
                    include_parent_locus=include_parent_locus
                )
                # If loci randomly overlap, resample
                while len(seen.intersection(candidates)) > 0:
                    candidates = self.bootstrap_candidate_loci(
                        locus, flank_limit=flank_limit,
                        window_size=window_size, chain=True,
                        include_parent_locus=include_parent_locus
                    )
                # Add all new bootstrapped loci to the seen list
                seen |= set(candidates)
                bootstraps.append(candidates)
            if chain:
                bootstraps = list(seen)
            #self.log.info("Found {} bootstraps",len(bootstraps))
            return bootstraps


    def pairwise_distance(self, loci=None):
        '''
            returns a vector containing the pairwise distances between loci
            in loci in vector form. See np.squareform for matrix
            conversion.
        '''
        import pandas as pd
        if loci is None:
            loci = list(self.iter_loci())
        query = '''
                SELECT id, chromosome, start, end FROM loci
                WHERE id in ("{}")
                ORDER BY id
        '''.format('","'.join([x.id for x in loci]))
        # extract chromosome row ids and start positions for each locus
        positions = pd.DataFrame(
            # Grab the chromosomes rowid because its numeric
            self._db.cursor().execute(query).fetchall(),
            columns=['locus','chrom','start','end']
        ).sort_values(by='locus')
        # chromosome needs to be floats
        positions.chrom = positions.chrom.astype('float')
        # Do a couple of checks
        assert len(positions) == len(loci), \
            'Some genes in dataset not if RefGen'
        assert all(positions.gene == [g.id for g in loci]), \
            'Genes are not in the correct order!'
        distances = LocusDist.gene_distances(
            positions.chrom.values,
            positions.start.values,
            positions.end.values
        )
        return distances

    def summary(self):
        print ("\n".join([
            'Loci: {} ',
            '{} sites']).format(
                self.name,self.num_loci(),
            )
        )

    '''----------------------------------------------------------------
        Aliases
    '''

    def add_alias(self, locus_id, alias):
        '''
            Add an alias name for a locus

            Parameters
            ----------
            locus_id: str
                The id of the locus that has an alias
            alias: str
                The alias of the locus
        '''
        if locus_id not in self:
            raise ValueError(f'{locus_id} not in the database')
        cur = self._db.cursor()
        cur.execute(
            'INSERT INTO aliases (id,alias) VALUES (?,?)',
            (locus_id,alias)
        )
        return True

    def num_aliases(self):
        '''
            Returns the number of aliases currently in the database

            Parameters
            ----------
            None

            Returns
            -------
            int : the number of aliases in the database

        '''
        num = self._db.cursor() \
            .execute('SELECT COUNT(*) FROM aliases') \
            .fetchone()
        if num == None:
            return 0
        else:
            return num[0]

    def _remove_aliases(self):
        '''
            Running this will delete the aliases from the 
            database. Warning! This is deletorious!

            Parameters
            ----------
            None

            Returns
            -------
            None
        '''
        self._db.cursor().execute('DELETE FROM aliases;')

    def add_annotations(self, filename, sep="\t", locus_col=0, skip_cols=None):
        ''' 
            Imports Annotation relationships from a csv file. By default will
            assume locus names are first column

            Parameters
            ----------
            filename : str 
                The file containing the annotations
            sep : str (default: \\t)
                The delimiter for the columns in the annotation file
            locus_col : int (default: 0)
                The index of the column containing the locus IDs
            skip_cols : default:None
                Optional names of columns to drop before adding 
                annotations

            Returns
            -------
            None if successful

        '''
        # import from file, assume right now that in correct order
        tbl = pd.read_table(filename,sep=sep,dtype=object)
        idx_name = tbl.columns[locus_col]
        tbl[idx_name] = tbl[idx_name].str.upper()
        # Set thie index to be the specified locus column
        tbl.set_index(idx_name,inplace=True)        
        
        # Drop columns if we need to
        if skip_cols is not None:
            # removing certain columns
            tbl.drop(tbl.columns[skip_cols],axis=1,inplace=True)
        
        # Get rid of any locus not in the refence genome
        cur = self._db.cursor()
        cur.execute('SELECT id FROM loci;')
        rm = set(tbl.index.values) - set([id[0] for id in cur.fetchall()])
        tbl.drop(rm,axis=0,inplace=True)
        del rm, cur
        
        # One Annotation per row, drop the nulls and duplicates
        tbl = tbl.reset_index()
        tbl = pd.melt(tbl,id_vars=idx_name,var_name='col',value_name='desc')
        tbl.drop('col',axis=1,inplace=True)
        tbl.dropna(axis=0,inplace=True)
        tbl.drop_duplicates(inplace=True)
        
        # Run the transaction to throw them in there
        cur = self._db.cursor()
        try:
            cur.execute('BEGIN TRANSACTION')
            cur.executemany(
                'INSERT INTO func VALUES (?,?)'
                ,tbl.itertuples(index=False))
            cur.execute('END TRANSACTION')
        
        except Exception as e:
            self.log.info("import failed: {}",e)
            cur.execute('ROLLBACK')
        
        # Make sure the indices are built
        self._build_indices()

    def remove_annotations(self):
        self._db.cursor().execute('DELETE FROM func;')

    def has_annotations(self) :
        cur = self._db.cursor()
        cur.execute('SELECT count(*) FROM func;')
        return (int(cur.fetchone()[0]) > 0)

    def get_annotations(self, item):
        # Build the query from all the loci provided
        if isinstance(item,(set,list)):
            ls = "{}".format("','".join([str(x) for x in item]))
            single = False
        else:
            ls = item
            single = True
        query = "SELECT * FROM func WHERE id IN ('{}');".format(ls)
        
        # Run the query and turn the result into a list of tuples
        cur = self._db.cursor()
        cur.execute(query)
        annotes = cur.fetchall()
        
        # If a list of loci was passed in, return a dictionary of lists
        if not single:
            res = {}
            for id,desc in annotes:
                if id in res:
                    res[id].append(desc)
                else:
                    res[id] = [desc]
        
        # Otherwise just return the list annotations
        else:
            res = []
            for id,desc in annotes:
                res.append(desc)
        return res
    def export_annotations(self, filename=None, sep="\t"):
        '''
            Make a table of all functional annotations.
        '''
        # Find the default filename
        if filename == None:
            filename = self.name + '_func.tsv'
        
        # Pull them all from sqlite
        cur = self._db.cursor()
        cur.execute("SELECT * FROM func;")
        
        # Used pandas to save it
        df = pd.DataFrame(cur.fetchall(),columns=['locus','desc']).set_index('locus')
        df.to_csv(filename,sep=sep)

    @classmethod
    def from_DataFrame(cls,df,name,description,build,organism,
            chrom_col='chrom',start_col='start',stop_col='stop',
            id_col='ID'):
        '''
            Imports a RefGen object from a CSV.
        '''
        self = cls.create(name,description,type='RefGen')
        self._global('build',build)
        self._global('organism','organism')
        loci = list()
        for i,row in df.iterrows():
            loci.append(
                Locus(
                    row[chrom_col],int(row[start_col]),int(row[stop_col]),
                    id=row[id_col]
                ).update(dict(row.items()))
            )
        self.add_loci(loci)
        self._build_indices()
        return self

    def _initialize_tables(self):
        '''
            Initializes the Tables holding all the information
            about the Loci. 
        '''
        cur = self._db.cursor()
        cur.execute('''
            CREATE TABLE IF NOT EXISTS loci (
                id TEXT NOT NULL UNIQUE,
                chromosome TEXT NOT NULL,
                start INTEGER,
                end INTEGER
            );
            /*
            ----Create a table that contains loci attribute
                mappings
            */
            CREATE TABLE IF NOT EXISTS loci_attrs (
                id TEXT NOT NULL,
                key TEXT,
                val TEXT,
                FOREIGN KEY(id) REFERENCES loci(id)
            );
            CREATE TABLE IF NOT EXISTS aliases (
                alias TEXT UNIQUE,
                id TEXT,
                FOREIGN KEY(id) REFERENCES loci(id)
            );
            CREATE TABLE IF NOT EXISTS func (
                id TEXT,
                desc TEXT,
                UNIQUE(id,desc) ON CONFLICT IGNORE,
                FOREIGN KEY(id) REFERENCES loci(id)
            );
            CREATE TABLE IF NOT EXISTS ortho_func (
                id TEXT,
                desc TEXT,
                UNIQUE(id,desc) ON CONFLICT IGNORE,
                FOREIGN KEY(id) REFERENCES loci(id)
            );
            
            CREATE INDEX IF NOT EXISTS loci_start_end ON loci (chromosome,start DESC, end ASC, id);
            CREATE INDEX IF NOT EXISTS loci_end_start ON loci (chromosome,end DESC,start DESC,id);
            CREATE INDEX IF NOT EXISTS loci_start ON loci (chromosome,start);
            CREATE INDEX IF NOT EXISTS loci_end ON loci (chromosome,end);
            CREATE INDEX IF NOT EXISTS locusid ON loci (id);
            CREATE INDEX IF NOT EXISTS locusattr ON loci_attrs (id);
            CREATE INDEX IF NOT EXISTS id ON func(id);
            CREATE INDEX IF NOT EXISTS id ON ortho_func(id); 
            ''');
