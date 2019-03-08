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
from functools import lru_cache

from .LocusDist import LocusDist
from .Locus import Locus
from .Exceptions import ZeroWindowError


class RefLoci(Freezable):
    """
        RefLoci are more than the sum of their parts. They have a name and
        represent something bigger than theirselves. They are important. They
        live on the disk in a database.
    """

    # Create a class-wide logger
    log = logging.getLogger(__name__)
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(asctime)s %(name)-12s %(levelname)-8s %(message)s")
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(logging.INFO)

    # Methods
    def __init__(self, name, basedir=None):
        super().__init__(name, basedir=basedir)
        self.name = name
        self._initialize_tables()

    def __len__(self):
        """
            Returns the number of loci in the dataset
        """
        return self._db.cursor().execute(""" SELECT COUNT(*) FROM loci""").fetchone()[0]

    # ---- Updated ----

    def __contains__(self, obj):
        # TODO 
        """
            Flexible on what you pass into the 'in' function
        """
        if isinstance(obj, Locus):
            # you can pass in a Locus object (this expression
            # should ALWAYS be true if you
            # created Locus object from this RefLoci object)
            if (
                self._db.cursor()
                .execute("""SELECT COUNT(*) FROM loci WHERE id = ?""", (obj.id,))
                .fetchone()[0]
                == 1
            ):
                return True
            else:
                return False
        elif isinstance(obj, str):
            # Can be a string object
            if (
                self._db.cursor()
                .execute(
                    """
                SELECT COUNT(*) FROM loci WHERE id = ?""",
                    (obj,),
                )
                .fetchone()[0]
                == 1
            ):
                return True
            elif (
                self._db.cursor()
                .execute(
                    """
                SELECT COUNT(*) FROM aliases WHERE alias = ?""",
                    (obj,),
                )
                .fetchone()[0]
                == 1
            ):
                return True
            else:
                return False
        elif isinstance(obj, tuple):
            chrom, start, end = obj
            if (
                self._db.cursor()
                .execute(
                    """
                SELECT COUNT(*) FROM loci WHERE chromosome=? AND start=? and end=?
                """,
                    (chrom, start, end),
                )
                .fetchone()[0]
                > 0
            ):
                return True
            else:
                return False
        else:
            raise TypeError("Cannot test for containment for {}".format(obj))

    def __getitem__(self, item):
        # TODO 
        """
            A convenience method to extract loci from the reference genome.
        """
        # NOTE: Dont LRU cache this, it gets cached in from_id
        if isinstance(item, str):
            return self.get_locus_from_id(item)
        elif isinstance(item, slice):
            # Map the slice vocab to the locus vocab
            chrom = item.start
            start = item.stop
            end = item.step
            return self.loci_within(Locus(chrom, start, end))
        # Allow for iterables of items to be passed in
        elif isinstance(item, tuple):
            chrom, start, end = tuple
            return self.loci_within(Locus(chrom, start, end))

    @lru_cache(maxsize=131072)
    def get_locus_from_id(self, locus_id):
        # TODO 
        """
            Returns a locus object from a string

            Parameters
            ----------
            locus_id : str
                ID of the locus you want to pull out

            Returns
            -------
            A single Locus Object

        """
        cur = self._db.cursor()
        # The most common use case is a direct reference to an id
        info = cur.execute(
            """
                SELECT chromosome,start,end,id,feature_type,score,strand,frame,source 
                FROM loci 
                WHERE id = ?
                """,
            [locus_id],
        ).fetchone()
        if info is None:
            # Try to fetch an alias
            (locus_id,) = (
                self._db.cursor()
                .execute(
                    """
                    SELECT id FROM aliases
                    WHERE alias = ?
                    """,
                    (locus_id,),
                )
                .fetchone()
            )
            return self.get_locus_from_id(locus_id)
        else:
            locus = Locus(*info)
            # Fetch the attrs
            attrs = cur.execute(
                "SELECT key,val FROM loci_attrs WHERE id = ?;", (locus_id,)
            )
            for key, val in attrs:
                locus[key] = val
            return locus


    def __iter__(self):
        # TODO 
        return self.iter_loci()
    
    #-------------------------------------------------------
    #            Putters
    #-------------------------------------------------------

    def add_locus(self, locus):
        import warnings
        warnings.warn(
            "The add_locus method is deprecated, a single locus can now be added "
            "using the add_loci method"
        )
        self.add_loci(locus)

    def add_loci(self, loci):
        # TODO 
        """
            Add loci to the database. Will also handle a single locus by wrapping it 
            in a list for you.

            Parameters
            ----------
            loci : an iterable of loci
                These loci are added to the database

            Returns
            -------
            None
        """
        if isinstance(loci, Locus):
            loci = [loci]
        try:
            # support adding lists of loci
            self.log.info("Adding {} RefLoci info to database".format(len(loci)))
            cur = self._db.cursor()
            cur.execute("BEGIN TRANSACTION")
            cur.executemany(
                """
                INSERT INTO loci 
                (id,chromosome,start,end,feature_type,score,strand,frame,source) 
                VALUES (?,?,?,?,?,?,?,?,?)
                """,
                ((x.name, x.chrom, x.start, x.end, x.feature_type, x.score, x.strand, x.frame, x.source)\
                 for x in loci),
            )
            breakpoint()
            cur.executemany(
                "INSERT OR REPLACE INTO loci_attrs (id,key,val) VALUES (?,?,?)",
                ((x.id, key, val) for x in loci for key, val in x.attr.items()),
            )
            cur.execute("END TRANSACTION")
            self._update_cache()
        except Exception as e:
            cur.execute("ROLLBACK")
            raise e

    def remove_locus(self, item):
        # TODO 
        """
            Remove loci from the database

            Parameters
            ----------
            loci : an iterable of IDs or a single ID

            Returns
            -------
            The number of loci removed
        """
        if isinstance(item, str):
            # We have a single locus
            cur = self._db.cursor()
            cur.execute("DELETE FROM loci WHERE id = ?", (item,))
        elif isinstance(item, tuple):
            # We have coordinates
            pass
        else:
            raise ValueError("Cannot find: {}. Must be an ID or coordinate tuple.")

    def import_gff(self, filename, feature_type="*", ID_attr="Name", attr_split="="):
        # TODO 
        """
            Imports RefLoci from a gff (General Feature Format) file.
            See more about the format here:
            http://www.ensembl.org/info/website/upload/gff.html

            Parameters
            ----------

            filename : str
                The path to the GFF file.
            name : str
                The name if the RefLoci object to be stored in the core
                minus80 database.
            feature_type : str (default: '*')
                The name of the feature (in column 2) that designates a
                locus. These features will be the main object that the RefLoci
                encompasses. The default value '*' will import all feature
                types.
            ID_attr : str (default: ID)
                The key in the attribute column which designates the ID or
                name of the feature.
            attr_split : str (default: '=')
                The delimiter for keys and values in the attribute column
        """
        loci = list()
        chroms = dict()
        if filename.endswith(".gz"):
            IN = gzip.open(filename, "rt")
        else:
            IN = open(filename, "r")
        for line in IN:
            # skip comment lines
            if line.startswith("#"):
                continue
            # Get the main information
            (
                chrom,
                source,
                feature,
                start,
                end,
                score,
                strand,
                frame,
                attributes,
            ) = line.strip().split("\t")
            # Get the attributes
            attributes = dict(
                [
                    (field.strip().split(attr_split))
                    for field in attributes.strip(";").split(";")
                ]
            )
            # If feature_type option is '*' short-circuit and import it
            if feature_type == '*' or feature == feature_type:
                loci.append(
                    Locus(
                        chrom,
                        int(start),
                        int(end),
                        attributes[ID_attr].upper().strip('"'),
                        strand=strand,
                        **attributes,
                    ).update(attributes)
                )
                # Keep track of seen chromosomes
                if chrom not in chroms:
                    chroms[chrom] = end
                else:
                    if end > chroms[chrom]:
                        chroms[chrom] = end
        IN.close()
        return loci
        self.add_loci(loci)

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

    def iter_loci(self):
        # TODO 
        """
            Iterates over loci in database.

            Returns
            -------
            A generator containing loci
        """
        for (id,) in self._db.cursor().execute(
            """
                SELECT id FROM loci
            """
        ):
            yield self[id]

    def random_locus(self, **kwargs):
        # TODO 
        """
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
            A Locus object (locuspocus.Locus based)

        """
        (id,) = (
            self._db.cursor()
            .execute(
                """
            SELECT id from loci ORDER BY RANDOM() LIMIT 1;
        """
            )
            .fetchone()
        )
        return self.get_locus_from_id(id, **kwargs)

    def random_loci(self, n, **kwargs):
        # TODO 
        """
            Return random loci from the database, without replacement.

            Parameters
            ----------
            n : int

            **kwargs : key,value pairs
                Extra parameters passed onto the locus init method

            Returns
            -------
            An iterable containing n (unique) random loci

        """
        return [
            self.get_locus_from_id(x, **kwargs)
            for (x,) in self._db.cursor().execute(
                "SELECT id from loci ORDER BY RANDOM() LIMIT ?", (n,)
            )
        ]

    def intersection(self, loci):
        # TODO 
        """
            Return the subset of loci that are in the dataset.

            Parameters
            ----------
            loci : list-like of Locus

            Returns
            -------
            a list like object containing loci that
            are in the dataset.
        """
        return [x for x in loci if x in self]


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

    """----------------------------------------------------------------
        Aliases
    """

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

    def _remove_aliases(self):
        # TODO 
        """
            Running this will delete the aliases from the
            database. Warning! This is deletorious!

            Parameters
            ----------
            None

            Returns
            -------
            None
        """
        self._db.cursor().execute("DELETE FROM aliases;")

    def _update_cache(self):
        # TODO 
        self.num_loci.cache_clear()

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
                source TEXT,
                feature_type TEXT,
                start INTEGER,
                end INTEGER,

                /* Add in the rest of the GFF fields  */
                score FLOAT,
                strand INT,
                frame INT
                
                /* Add in a constraint  */
                UNIQUE(chromosome,source,feature_type,start,end) 
            );
            '''
        )
        cur.execute(
        # Create a table that contains loci attribute mapping
        '''
            CREATE TABLE IF NOT EXISTS loci_attrs (
                LID INT NOT NULL,
                key TEXT,
                val TEXT,
                FOREIGN KEY(LID) REFERENCES loci(LID)
            );
            '''
        )
        cur.execute(
        # Create a table with aliases
        '''
            CREATE TABLE IF NOT EXISTS aliases (
                alias TEXT PRIMARY KEY,
                LID INTEGER,
                FOREIGN KEY(LID) REFERENCES loci(LID)
            );
        ''')
        cur.execute(
        # Create a view with names
            '''
            CREATE VIEW IF NOT EXISTS named_loci AS
            SELECT alias, chromosome, source, feature_type, start, end, score, strand, frame
            FROM aliases 
            JOIN loci ON aliases.LID = loci.LID;
            '''
        )
        cur.execute(
            '''
            CREATE TRIGGER IF NOT EXISTS assign_LID INSTEAD OF INSERT ON named_loci
            FOR EACH ROW
            BEGIN
                INSERT OR IGNORE INTO loci 
                (chromosome,source,feature_type,start,end,score,strand,frame)
                VALUES 
                (NEW.chromosome, NEW.source, NEW.feature_type, NEW.start, NEW.end, NEW.score, NEW.strand, NEW.frame);
                INSERT INTO aliases 
                SELECT NEW.alias, last_insert_rowid();
            END
            '''
        )
