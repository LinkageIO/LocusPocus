#!/usr/bin/python3
import itertools
import random
import scipy as sp
import gzip
import logging

import numpy as np

from minus80 import Freezable
from collections.abc import Iterable

from .Locus import Locus
from .Exceptions import ZeroWindowError


class RefLoci(Freezable):
    """
        RefLoci are more than the sum of their parts. They have a name and
        represent something bigger than theirselves. They are important. They
        live on the disk in a database.
    """

    # Create a class-wide logger

    # Methods
    def __init__(self, name, basedir=None):
        super().__init__(name, basedir=basedir)
        self.name = name
        self._initialize_tables()

        self.log = logging.getLogger(__name__)
        handler = logging.StreamHandler()
        formatter = logging.Formatter("%(asctime)s %(name)-12s %(levelname)-8s %(message)s")
        handler.setFormatter(formatter)
        self.log.addHandler(handler)
        self.log.setLevel(logging.INFO)


    def __len__(self):
        """
            Returns the number of loci in the dataset
        """
        return self._db.cursor().execute(
            """ SELECT COUNT(*) FROM loci"""
        ).fetchone()[0]

    def _get_locus_by_LID(self,LID):
        '''
            Get a locus by its LID
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
                refloci=self
            )
            return locus
        except TypeError as e:
            raise ValueError(f'There is no locus in the database with LID of {LID}')

    def _get_LID(self,locus,cur=None):
        '''
            Return the Locus Identifier used internally by RefLoci

            Parameters
            ----------
            locus : one of (str,Locus)
                The locus for which to find the LID of, this can
                be either a Locus object OR a name/alias
            cur : apsw Cursor object
                If you are looking up LIDS within a transaction
                you need to pass in the cursor you are using or
                as the LIDS are not visible to new cursors until
                the transaction has completed.
        '''
        if cur is None:
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

    def add_locus(self, locus, cur=None):
        """
            Add locus to the database. 

            Parameters
            ----------
            locus : a Locus object
                This locus will be added to the db
            cur : a db curson
                An optional cursor object to use. If none, a 
                cursor will be created. If importing many 
                loci in a loop, use a bulk transaction.

            Returns
            -------
            The locus ID of the feshly added locus
        """
        if cur is None:
            cur = self._db.cursor()

        if locus.name == 'AC148152.3_FGT001':
            print(f"The hash for {locus.name} is {hash(locus)}")
            print(f"{locus}")
        # insert the core feature data
        core,attrs = locus.as_record()
        cur.execute(
            """
            INSERT OR IGNORE INTO loci 
                (chromosome,start,end,source,feature_type,strand,frame,name,hash)
                VALUES (?,?,?,?,?,?,?,?,?)
            """,
            core,
        )
        # get the fresh LID
        (LID,) = cur.execute('SELECT last_insert_rowid()').fetchone()
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
        filename, 
        feature_type="*", 
        ID_attr="ID", 
        parent_attr='Parent',
        attr_split="="
    ):
        """
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
        """
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
                del attributes[ID_attr]
            else:
                name = None
            # Parse out the parent info
            if parent_attr in attributes:
                parent = attributes[parent_attr]
                del attributes[parent_attr]
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
        self.log.info(f'Found {len(loci)} primary loci and {total_loci} total_loci, adding them to database')
        IN.close()
        with self._bulk_transaction() as cur:
            for l in loci:
                self.add_locus(l,cur=cur)
        self.log.info('Done!')

    def __contains__(self, locus):
        """
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
        """
        try:
            # If we can get an LID, it exists
            LID = self._get_LID(locus)
            return True
        except ValueError as e:
            return False

    def __getitem__(self, item):
        """
            A convenience method to extract loci from the reference genome.
        """
        LID = self._get_LID(item)
        return self._get_locus_by_LID(LID)

    def __iter__(self):
        raise NotImplementedError()
    
    def __delitem__(self,item):
        """
            Remove a locus and all of its sub loci from the database

            Parameters
            ----------
            loci : an iterable of IDs or a single ID

            Return: bool (default:True)
                By default, the loci need to be distinct so 
                no duplicates are allowed in the output. If 
                set to False, this check will not be made 
                and duplicates are possible

            Returns
            -------
            A list of n Locus objects

        """
        loci = set()
        maxLID, = self._db.cursor().execute(
            'SELECT MAX(LID) from loci;'
        ).fetchone()
        if n > maxLID:
            raise ValueError('More than the maximum loci in the database was requested')
        while len(loci) < n:
            try:
                rand_LID = int(np.random.randint(0,maxLID,1)[0])
                loci.add(self._get_locus_by_LID(rand_LID))
            except ValueError as e:
                #self.log.info(f'hit a tombstone! LID {rand_LID}')
                continue
        return list(loci)


    # ----------- Internal Methods  -----------#


    def summarize_feature_types(self):
        # TODO 
        """
            Helper to summarize the Loci by feature type
        """
        return (
            self._db.cursor()
            .execute(
                """ 
           SELECT feature_type, COUNT(feature_type) 
           FROM loci 
           GROUP BY feature_type 
           ORDER BY COUNT(feature_type) DESC
           """
            )
            .fetchall()
        )


    def upstream_loci(
        # TODO 
        self, locus, locus_limit=1000, window_size=None, within=False, feature_type="%"
    ):
        """
            Find loci that START upstream of a locus.
            Loci are ordered so that the nearest loci are
            at the beginning of the list.

            Return Loci that overlap with the upstream window,
            This includes partially overlapping loci, but NOT
            loci that are returned by the loci_within method.

            Looks like: (y=yes,returned; n=no,not returned)

            nnn  yyyyyyy   yyyyyy   yyyyy  yyyyyy nnnn nnnn nnnnnnnn
                                             nnnn
                                           start             end
                -----------------------------x****************x--
                   ^_________________________| Window (upstream)
        """
        if locus.window == 0 and window_size is None:
            raise ValueError(
                "Asking for upstream loci for {} which has a window size of 0".format(
                    locus
                )
            )
        if window_size is not None:
            upstream = locus.start - window_size
        else:
            upstream = locus.upstream
        if within:
            return [
                self.get_locus_from_id(x)
                for (x,) in self._db.cursor().execute(
                    """
                    SELECT id FROM loci
                    INDEXED BY loci_start_end
                    WHERE chromosome = ?
                    AND start >= ?  -- Gene must end AFTER locus window (upstream)
                    AND end < ? -- Gene must end BEFORE locus starts
                    AND feature_type LIKE ?
                    ORDER BY start DESC
                    LIMIT ?
                    """,
                    (locus.chrom, upstream, locus.start, feature_type, locus_limit),
                )
            ]
        else:
            return [
                self.get_locus_from_id(x)
                for (x,) in self._db.cursor().execute(
                    """
                    SELECT id FROM loci
                    INDEXED BY loci_start_end
                    WHERE chromosome = ?
                    AND end >= ?  -- Gene must end AFTER locus window (upstream)
                    AND end < ? -- Gene must end BEFORE locus ends
                    AND start < ? -- Gene must start before locus starts
                    AND feature_type LIKE ?
                    ORDER BY start DESC
                    LIMIT ?
                    """,
                    (
                        locus.chrom,
                        upstream,
                        locus.end,
                        locus.start,
                        feature_type,
                        locus_limit,
                    ),
                )
            ]

    def downstream_loci(
        # TODO 
        self, locus, locus_limit=1000, window_size=None, within=False, feature_type="%"
    ):
        """
            Returns loci downstream of a locus. Loci are ordered
            so that the nearest loci are at the beginning of the list.

            Return Loci that overlap with the downstream window,
            This includes partially overlapping loci, but NOT
            loci that are returned by the loci_within method.

            Looks like: (y=yes,returned; n=no,not returned)

            nnn  nnnnnnn   nnnnnn nnnn  yyyy  yyyyyy yyyy yyyyyy  nnnnn
               start             end
              ---x****************x--------------------------------
                                  |_______________________^ Window (downstream)
        """
        if locus.window == 0 and window_size is None:
            raise ValueError(
                "Asking for upstream loci for {} with no window size", locus.id
            )
        if window_size is not None:
            downstream = locus.end + window_size
        else:
            downstream = locus.downstream
        if within:
            return [
                self.get_locus_from_id(x)
                for (x,) in self._db.cursor().execute(
                    """
                    SELECT id FROM loci
                    INDEXED BY loci_start_end
                    WHERE chromosome = ?
                    AND start >= ?
                    AND end <= ?
                    AND feature_type LIKE ?
                    ORDER BY start ASC
                    LIMIT ?
                    """,
                    (locus.chrom, locus.end, downstream, feature_type, locus_limit),
                )
            ]
        else:
            return [
                self.get_locus_from_id(x)
                for (x,) in self._db.cursor().execute(
                    """
                    SELECT id FROM loci
                    INDEXED BY loci_start_end
                    WHERE chromosome = ?
                    AND start > ?
                    AND start <= ?
                    AND end > ?
                    AND feature_type LIKE ?
                    ORDER BY start ASC
                    LIMIT ?
                    """,
                    (
                        locus.chrom,
                        locus.start,
                        downstream,
                        locus.end,
                        feature_type,
                        locus_limit,
                    ),
                )
            ]

    def flanking_loci(
        # TODO 
        self,
        loci,
        flank_limit=1,
        chain=True,
        window_size=None,
        within=False,
        feature_type="%",
    ):
        """
            Returns loci upstream and downstream from a locus
            ** done NOT include loci within locus **
        """
        if isinstance(loci, Locus):
            # If we cant iterate, we have a single locus
            locus = loci
            if locus.window == 0 and window_size is None:
                raise ValueError(
                    "Asking for upstream loci for {} and no window size.", locus.id
                )
            upstream_locus_limit = int(flank_limit)
            downstream_locus_limit = int(flank_limit)
            up_loci = self.upstream_loci(
                locus,
                locus_limit=upstream_locus_limit,
                window_size=window_size,
                within=within,
                feature_type=feature_type,
            )
            down_loci = self.downstream_loci(
                locus,
                locus_limit=downstream_locus_limit,
                window_size=window_size,
                within=within,
                feature_type=feature_type,
            )
            if chain:
                return list(itertools.chain(up_loci, down_loci))
            return (up_loci, down_loci)
        else:
            iterator = iter(loci)
            loci = [
                self.flanking_loci(
                    locus,
                    flank_limit=flank_limit,
                    window_size=window_size,
                    within=within,
                    feature_type=feature_type,
                )
                for locus in iterator
            ]
            if chain:
                loci = list(itertools.chain(*loci))
            return loci

    def captured(self, loci, chain=True, feature_type="%"):
        # TODO 
        """
            Returns the Loci captured by the locus. This method returns
            loci that are *compeltely* within the specified loci.
            It will NOT return loci that are partially overlapping that
            are return using the within method.

            Parameters
            ----------
            loci : an iterable of loci
                The method will return encompassing loci for each
                locus in ther iterable. If a single locus is passed,
                a single result will be returned.
            chain : bool (default: True)
                Specifies whether or not to chain results. If true
                a single list will be returned, if False, a result for
                each locus passed in will be returned.
            feature : str (default: %)
                Specifies whether any locus in the database meeting the
                requirements will be returned or only loci of a certain
                type (gene, transposable_element, snp, etc.)
        """
        if isinstance(loci, Locus):
            return [
                self.get_locus_from_id(x)
                for (x,) in self._db.cursor().execute(
                    """
                    SELECT id FROM loci
                    WHERE chromosome = ?
                    AND start >= ?
                    AND end <= ?
                    AND feature_type LIKE ?
                    """,
                    (loci.chrom, loci.start, loci.end, feature),
                )
            ]
        else:
            iterator = iter(loci)
            loci = [self.loci_within(locus, chain=chain) for locus in iterator]
            if chain:
                loci = list(itertools.chain(*loci))
            return loci

    def encompassing_loci(self, loci, chain=True):
        # TODO 
        """
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
        """
        if isinstance(loci, Locus):
            return [
                self.get_locus_from_id(x)
                for (x,) in self._db.cursor().execute(
                    """
                    SELECT id FROM loci
                    WHERE chromosome = ?
                    AND start <= ? AND end >= ?
                    """,
                    (loci.chrom, loci.start, loci.end),
                )
            ]
        else:
            iterator = iter(loci)
            loci = [self.encompassing_loci(locus, chain=chain) for locus in iterator]
            if chain:
                loci = list(itertools.chain(*loci))
            return loci

    def loci_within(self, loci, chain=True, feature="%"):
        # TODO 
        """
            Returns the loci that START or END within a locus
            start/end boundary.

            Looks like: (y=yes,returned; n=no,not returned)

            nnn  nnnnnnn  yyyyyy   yyyyy  yyyyyy yyyyyyy
                    start                        end
                -----x****************************x-------------

        """
        if isinstance(loci, Locus):
            return [
                self.get_locus_from_id(x)
                for (x,) in self._db.cursor().execute(
                    """
                    SELECT id FROM loci
                    WHERE chromosome = ?
                    AND start >= ? AND start <= ?
                    """,
                    (loci.chrom, loci.start, loci.end),
                )
            ]
        else:
            iterator = iter(loci)
            loci = [self.loci_within(locus, chain=chain) for locus in iterator]
            if chain:
                loci = list(itertools.chain(*loci))
            return loci

    def candidate_loci(
        # TODO 
        self,
        loci,
        flank_limit=2,
        chain=True,
        window_size=None,
        include_parent_locus=False,
        include_parent_attrs=False,
        include_num_intervening=False,
        include_rank_intervening=False,
        include_num_siblings=False,
        include_SNP_distance=False,
        attrs=None,
        return_table=False,
    ):
        """
            Locus to locus mapping.
            Return loci between locus start and stop, plus additional
            flanking loci (up to flank_limit)

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

        """
        if isinstance(loci, Locus):
            # If not an iterator, its a single locus
            locus = loci
            loci_within = self.loci_within(locus)
            up_loci, down_loci = self.flanking_loci(
                locus, flank_limit=flank_limit, chain=False, window_size=window_size
            )

            # This always returns candidates together, if
            # you want specific up,within and down loci
            # use the specific methods
            candidates = sorted(itertools.chain(up_loci, loci_within, down_loci))
            # include the number of effective loci
            if include_rank_intervening == True:
                ranks = sp.stats.rankdata(
                    [abs(x.center_distance(locus)) for x in candidates]
                )
            # Iterate through candidate loci and propagate the
            # parental info
            for i, cand in enumerate(candidates):
                # include parent locus id if thats specified
                if include_parent_locus == True:
                    cand.update({"parent_locus": locus.id})
                if include_rank_intervening == True:
                    cand.update({"intervening_rank": ranks[i]})
                # update all the parent_attrs
                if include_parent_attrs and len(include_parent_attrs) > 0:
                    if "all" in include_parent_attrs:
                        include_parent_attrs = locus.attr.keys()
                    for attr in include_parent_attrs:
                        attr_name = "parent_{}".format(attr)
                        cand.update({attr_name: locus[attr]})
            if include_num_intervening == True:
                num_down = 0
                num_up = 0
                # Sort the loci by their distance from the locus
                cands_with_distances = [
                    (cand, abs(cand.center_distance(locus))) for cand in candidates
                ]
                cands_with_distances = sorted(cands_with_distances, key=lambda x: x[1])
                for cand, distance in cands_with_distances:
                    if locus.within(cand):
                        cand.update({"num_intervening": -1})
                    elif cand.center >= locus.center:
                        cand.update({"num_intervening": num_down})
                        num_down += 1
                    elif cand.center <= locus.center:
                        cand.update({"num_intervening": num_up})
                        num_up += 1
            if include_num_siblings == True:
                for cand in candidates:
                    cand.update({"num_siblings": len(candidates)})
            if include_SNP_distance == True:
                for cand in candidates:
                    distance = abs(cand.center_distance(locus))
                    cand.update({"SNP_distance": distance})
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
                    attrs=attrs,
                )
                for locus in iterator
            ]
            if chain:
                if return_table:
                    candidates = pd.concat(candidates)
                else:
                    candidates = list(set(itertools.chain(*candidates)))
            return candidates

    def bootstrap_candidate_loci(
        # TODO 
        self,
        loci,
        flank_limit=2,
        chain=True,
        window_size=None,
        include_parent_locus=False,
    ):
        """
            Returns candidate loci which are random, but conserves
            total number of overall loci.

            Parameters
            ----------
            loci : locuspocus.Locus (also handles an iterable containing Loci)
                a locus or iterable of loci
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

        """
        if isinstance(loci, Locus):
            # We now have a single locus
            locus = loci
            # grab the actual candidate loci
            num_candidates = len(
                self.candidate_loci(
                    locus, flank_limit=flank_limit, chain=True, window_size=window_size
                )
            )
            if num_candidates == 0:
                return []
            # Snps a random loci from the genome
            random_locus = self.random_locus()
            # Snag the same number of candidates
            random_candidates = self.upstream_loci(
                random_locus, locus_limit=num_candidates, window_size=10e100
            )
            if len(random_candidates) != num_candidates:
                # somehow we hit the end of a chromosome
                # or something, just recurse
                random_candidates = self.bootstrap_candidate_loci(
                    locus, flank_limit=flank_limit, chain=True
                )
            if include_parent_locus == True:
                for cand in random_candidates:
                    cand.update({"parent_locus": random_locus.id})
            return random_candidates
        else:
            # Sort the loci so we can collapse down
            locus_list = sorted(loci)
            seen = set()
            bootstraps = list()
            target = self.candidate_loci(
                locus_list,
                flank_limit=flank_limit,
                chain=False,
                window_size=window_size,
            )
            target_accumulator = 0
            candidate_accumulator = 0
            # self.log.info('target: {}, loci: {}',len(target),len(locus_list))
            for i, (locus, targ) in enumerate(zip(locus_list, target)):
                # compare downstream of last locus to current locus
                candidates = self.bootstrap_candidate_loci(
                    locus,
                    flank_limit=flank_limit,
                    chain=True,
                    window_size=window_size,
                    include_parent_locus=include_parent_locus,
                )
                # If loci randomly overlap, resample
                while len(seen.intersection(candidates)) > 0:
                    candidates = self.bootstrap_candidate_loci(
                        locus,
                        flank_limit=flank_limit,
                        window_size=window_size,
                        chain=True,
                        include_parent_locus=include_parent_locus,
                    )
                # Add all new bootstrapped loci to the seen list
                seen |= set(candidates)
                bootstraps.append(candidates)
            if chain:
                bootstraps = list(seen)
            # self.log.info("Found {} bootstraps",len(bootstraps))
            return bootstraps


    def summary(self):
        # TODO 
        print("\n".join(["Loci: {} ", "{} sites"]).format(self.name, self.num_loci()))


    def add_alias(self, locus_id, alias):
        # TODO 
        """
            Add an alias name for a locus

            Parameters
            ----------
            locus_id: str
                The id of the locus that has an alias
            alias: str
                The alias of the locus
        """
        if locus_id not in self:
            raise ValueError(f"{locus_id} not in the database")
        cur = self._db.cursor()
        cur.execute("INSERT INTO aliases (id,alias) VALUES (?,?)", (locus_id, alias))
        return True

    def num_aliases(self):
        # TODO 
        """
            Returns the number of aliases currently in the database

            Parameters
            ----------
            None

            Returns
            -------
            int : the number of aliases in the database

        """
        num = self._db.cursor().execute("SELECT COUNT(*) FROM aliases").fetchone()
        if num == None:
            return 0
        else:
            return num[0]

    def get_feature_list(self, feature="%"):
        # TODO 
        '''
        Get a list of all of the loci in the database of the specified feature type
        '''
        cur = self._db.cursor()
        cur.execute(
            """
            SELECT id FROM loci
            WHERE feature_type LIKE ?
            ORDER BY chromosome, start
            """,
            (feature,),
        )
        p = cur.fetchall()
        feature_list = list(itertools.chain.from_iterable(p))
        return feature_list

    # -----------------------------------------
    #
    #           Class Methods
    #
    # -----------------------------------------
    
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
            '''
        )
        self._initialize_tables()

    def _initialize_tables(self):
        """
            Initializes the Tables holding all the information
            about the Loci.
        """
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
        # Create a view that inlcudes only top-level primary loci
        # In this table, the LIDS are loci that do not have any parents
        # and are thus top-level or "primary" features such
        cur.execute(
        '''
            CREATE TEMP VIEW IF NOT EXISTS primary_loci AS 
                SELECT LID FROM loci l 
                LEFT JOIN relationships r 
                ON l.LID = r.child 
                WHERE child IS NULL -- if the child listing is NULL, there is no parent locus
        '''
        )

