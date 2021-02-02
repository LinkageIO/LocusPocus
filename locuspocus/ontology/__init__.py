#!/usr/bin/python3

import copy
import logging
import re

import pandas as pd
import minus80 as m80

from minus80 import Freezable
from minus80.Tools import rawFile
from pandas import DataFrame
from scipy.stats import hypergeom
from functools import lru_cache
from collections import OrderedDict, defaultdict
from itertools import chain

from typing import Optional, Iterable, List

from .term import Term
from locuspocus import Loci
from locuspocus.exceptions import MissingLocusError

__all__ = ["Ontology", "Term"]

log = logging.getLogger(__name__)


class Ontology(Freezable):
    """
    An Ontology is just a collection of terms. Each term is just a
    collection of genes. Sometimes terms are related or nested
    within each other, sometimes not. Simple enough.

    Parameters
    ----------
    name : unique identifier

    Returns
    -------
    An Ontology Object

    """

    def __init__(self, name, rootdir: Optional[str] = None):
        super().__init__(name, rootdir=rootdir)
        self._initialize_tables()
        self.metadata = self.m80.doc.table("metadata")

        self._loci = None

    def __len__(self):
        """
        Return the number of non-empty terms
        """
        return self.num_terms(min_term_size=1)

    def __iter__(self):
        return self.terms()

    # -----------------------------------------
    #       Properties
    # -----------------------------------------

    @property
    def loci(self) -> Loci:
        # lazy evaluation
        if self._loci is None:
            self._loci = Loci(self.m80.name, rootdir=self.m80.thawed_dir)
        return self._loci

    # -----------------------------------------
    #       Methods
    # -----------------------------------------

    def add_term(self, term, cursor=None, loci_cursor=None):
        """
        This will add a single term to the ontology

        Parameters
        ----------
        term : Term object
            The term object you wish to add.
        cursor : apsw cursor object
            A initialized cursor object for the term db, for batch operation.
            This will allow for adding many terms in one transaction as long as
            the passed in cursor has executed the "BEGIN TRANSACTION" command.
        loci_cursor : apsw cursor object
            A initialized cursor object loci db, for batch operation. This will
            allow for adding many loci in one transaction as long as the
            passed in cursor has executed the "BEGIN TRANSACTION" command.
        """

        if not cursor:
            cur = self.m80.db.cursor()
            cur.execute("BEGIN TRANSACTION")
        else:
            cur = cursor

        if not loci_cursor:
            lcur = self.loci.m80.db.cursor()
            lcur.execute("BEGIN TRANSACTION")
        else:
            lcur = loci_cursor

        # Add the term id and description
        cur.execute(
            """
            INSERT OR ABORT INTO terms (name, desc)
            VALUES (?, ?)""",
            (term.name, term.desc),
        )

        (TID,) = cur.execute("SELECT last_insert_rowid()").fetchone()

        if TID is None:  # pragma: no cover
            # I dont know when this would happen without another exception being thrown
            raise ValueError(f"{term} was not assigned a valid TID!")

        for key, val in term.attrs.items():
            cur.executemany("""
                INSERT INTO term_attrs 
                    (TID, key, val) 
                    VALUES (?,?,?)
                """, ((TID,key,xval) for xval in val)
            )

        # separate the new loci from the existing loci
        new_LIDs = []
        existing_LIDs = []
        for l in term.loci:
            try:
                existing_LIDs.append(self.loci._get_LID(l))
            except MissingLocusError:
                new_LIDs.append(self.loci.add_locus(l, cur=lcur))

        for LID in new_LIDs + existing_LIDs:
            cur.execute(
                """
                INSERT INTO term_loci
                    (TID,LID)
                    VALUES (?,?)
            """,
                (TID, LID),
            )

        if not cursor:
            cur.execute("END TRANSACTION")

    def num_terms(self, min_term_size=0, max_term_size=10e10):
        """
        Returns the number of terms in the Ontology
        within the min_term_size and max_term_size

        Parameters
        ----------
        min_term_size (default:0)
            The minimum number of loci associated with the term
        max_term_size (default: 10e10)
            The maximum number of loci associated with the term

        Returns
        -------
        the number of terms that fit the criteria

        """
        return (
            self.m80.db.cursor()
            .execute(
                """SELECT COUNT(*) FROM (
                SELECT DISTINCT(TID) FROM term_loci 
                GROUP BY TID 
                HAVING COUNT(TID) >= ? 
                    AND  COUNT(TID) <= ?
            );""",
                (min_term_size, max_term_size),
            )
            .fetchone()[0]
        )

    @lru_cache(maxsize=65536)
    def __getitem__(self, item):
        """
        Retrieve a term by name or TID.
        """
        try:
            if isinstance(item, str):
                (TID, name, desc) = (
                    self.m80.db.cursor()
                    .execute(
                        "SELECT TID, name, desc from terms WHERE name = ?", (item,)
                    )
                    .fetchone()
                )
            elif isinstance(item, int):
                (TID, name, desc) = (
                    self.m80.db.cursor()
                    .execute("SELECT TID, name, desc from terms WHERE TID = ?", (item,))
                    .fetchone()
                )
            else:
                raise TypeError
        except TypeError as e:  # Not in database
            raise e
        else:
            term_loci = [
                self.loci._get_locus_by_LID(LID)
                for LID, in self.m80.db.cursor()
                .execute(""" SELECT LID FROM term_loci WHERE TID = ?""", (TID,))
                .fetchall()
            ]
            attrs = defaultdict(list)
            for k, v in self.m80.db.cursor().execute(
                """ SELECT key,val FROM term_attrs WHERE TId = ?""", (TID,)
            ):
                attrs[k].append(v)
            return Term(name=name, desc=desc, loci=term_loci, attrs=attrs)

    def terms_containing(
        self,
        loci: Loci,
        min_term_size: int = 0,
        max_term_size: int = 10e10,
    ):
        """
        Retrurns the set of terms which contains the
        specified loci.

        Parameters
        ----------
        loci : Loci
            The list of loci for which to retrieve
            corresponding terms.
        min_term_size : int (default: 0)
            The minimum term size for which to test enrichment. Useful
            for filtering out very small terms that would be uninformative
            (e.g. single gene terms)
        max_term_size : int (default: 10e10)
            The maximum term size for which to test enrichment. Useful
            for filtering out large terms that would otherwise be
            uninformative (e.g. top level GO terms)

        Returns
        -------
        list of terms which contain provided loci
        """
        # Extract LIDS for loci in Ontology
        LIDs = []
        for l in loci:
            try:
                LIDs.append(self.loci._get_LID(l))
            except MissingLocusError:
                continue
        # query the database
        TIDs = (
            self.m80.db.cursor()
            .execute("""
                SELECT DISTINCT TID 
                FROM term_loci WHERE LID IN ('{}')
            """.format("','".join(map(str, LIDs)))
            ).fetchall()
        )

        terms = [self[TID] for (TID,) in TIDs]
        # Fetch the terms with the proper size
        terms = list(
            filter(
                lambda t: (len(t) >= min_term_size) and (len(t) <= max_term_size), terms
            )
        )
        return terms

    def terms(self, min_term_size=0, max_term_size=10e10) -> Iterable[Term]:
        """
        Return a generator that iterates over each term in the ontology.
        """
        terms = self.m80.db.cursor().execute(
            """
            SELECT DISTINCT(TID) from term_loci
            GROUP BY TID
            HAVING COUNT(LID) >= ?
                AND COUNT(LID) <= ?
            """,
            (min_term_size, max_term_size),
        )
        for (id,) in terms:
            yield self[id]

    def summary(self) -> str:
        return "Ontology:{} -  contains {} terms containing {} distinct loci".format(
            self.m80.name, len(self), len(self.loci)
        )

    def rand(self, n=1, min_term_size=1, max_term_size=10e10, autopop=True) -> Term:
        """
        Return a random Term from the Ontology

        Parameters
        ----------
        n : int (default: 1)
            The number of random terms to return
        min_term_size : int (default: 1)
            The smallest acceptable term size
            i.e. the number of genes annotated to the term
        max_term_size : int (default: 10e10)
            The largest acceptable term size
        autopop: bool = True
            If a single term is found, return the term object,
            otherwise return a list of terms.
        """
        cur = self.m80.db.cursor()
        TIDs = cur.execute(
            """ 
            SELECT DISTINCT(TID) FROM term_loci 
            GROUP BY LID
            HAVING COUNT(LID) >= ?
                AND COUNT(LID) <= ?
            ORDER BY RANDOM() 
            LIMIT ?;
        """,
            (min_term_size, max_term_size, n),
        ).fetchall()
        if len(TIDs) == 0:
            raise ValueError(
                "No Terms exists with this criteria "
                "{} < len(term) < {}:".format(min_term_size, max_term_size)
            )
        terms = [self[TID] for TID, in TIDs]
        if len(terms) == 1 and autopop:
            return terms[0]
        else:
            return terms

    def parents(self, term, parent_attr='is_a'):
        """
            Return an iterable containing the parents of a term.
            Parents are determined via the is_a property of the term.

            Parameters
            ----------
            term : GOTerm

            Returns
            -------
            An iterable containing the parents
            NOTE: note guaranteed to be in order NOR guaranteed 
            to not contain duplicates
        """
        for parent in term.attrs[parent_attr]:
            yield self[parent]
            yield from self.parents(parent, parent_attr=parent_attr)

    def enrichment(
        self,
        loci,
        min_term_size=2,
        max_term_size=300,
        pval_cutoff=0.05,
        num_universe=None,
        return_table=False,
        label=None,
        include_genes=False,
        bonferroni_correction=True,
        min_overlap=1,
    ):
        """
        Evaluates enrichment of loci within the locus list for terms within
        the ontology. NOTE: this only tests terms that have at least one
        locus that exists in locus_list.

        Parameters
        ----------
        locus_list : list of Locus *of* instance of Ontology
            A list of loci for which to test enrichment. i.e. is there
            an over-representation of these loci within and the terms in
            the Ontology. If an ontology is passed, each term in the ontology
            will be iterated over and tested as if they were a locus_list.
        pval_cutoff : float (default: 0.05)
            Report terms with a pval lower than this value
        bonferroni_correction : bool (default: True)
            correct for testing multiple terms using Bonferroni correction
        max_term_size : int (default: 300)
            The maximum term size for which to test enrichment. Useful
            for filtering out large terms that would otherwise be
            uninformative (e.g. top level GO terms)
        min_term_size : int (default: 5)
            The minimum term size for which to test enrichment. Useful
            for filtering out very small terms that would be uninformative
            (e.g. single gene terms)
        num_universe : int (default: None)
            Use a custom universe size for the hypergeometric calculation,
            for instance if you have a reduced number of genes in a reference
            co-expression network. If None, the value will be calculated as
            the total number of distinct genes that are observed in the
            ontology.
        include_genes : bool (default: False)
            Include comma delimited genes as a field
        return_table : bool (default: False)
            If True, return results as a data frame
        label: str (default: None)
            If a label is specified, it will be inlcuded in the results
        min_overlap : int (default: 1)
            The minimum overlap between genes in the term and genes in
            the locus list. Increasing this value can minimize spurious
            or uninformative terms
        """
        if isinstance(loci, Ontology):
            ontology = loci
            log.info("Calculating enrichment for an  Ontology: {}", ontology.m80.name)
            enrich = []
            for term in ontology:
                e = self.enrichment(
                    term.loci,
                    pval_cutoff=pval_cutoff,
                    max_term_size=max_term_size,
                    min_term_size=min_term_size,
                    num_universe=num_universe,
                    return_table=return_table,
                    label=term.id,
                    include_genes=include_genes,
                    bonferroni_correction=bonferroni_correction,
                    min_overlap=min_overlap,
                )
                enrich.append(e)
            if return_table:
                return pd.concat(enrich)
            else:
                return enrich
        # return a new copy of each
        terms = [
            copy.copy(term)
            for term in self.terms_containing(
                loci, min_term_size=min_term_size, max_term_size=max_term_size
            )
        ]
        # Calculate the size of the Universe
        if num_universe is None:
            num_universe = self.num_distinct_loci()
        log.info(
            "{}: Loci occur in {} terms, containing {} genes".format(
                label, len(terms), num_universe
            )
        )
        significant_terms = []
        for term in terms:
            term_genes = set(term.loci)
            if len(term_genes) > max_term_size:
                continue
            num_common = len(term_genes.intersection(loci))
            num_in_term = len(term_genes)
            num_sampled = len(loci)
            # the reason this is num_common - 1 is because we are looking for 1 - cdf
            # and we need to greater than OR EQUAL TO num_common
            # Look. Do this in ipython:
            """
                In [99]: probs = [hypergeom.pmf(x,100,5,10) for x in range(0,6)]
                In [100]: probs
                Out[100]: 
                [0.58375236692612187,
                 0.33939091100357333,
                 0.070218809173150043,
                 0.006383528106649855,
                 0.00025103762217164457,
                 3.3471682956218215e-06]
                In [103]: 1-sum(probs[0:3]) 
                # Get the probs of drawing 3 or more
                Out[103]: 0.006637912897154763
                # Remember slicing is exclusive for the end value
                In [105]: hypergeom.sf(3,100,5,10)
                # That aint right
                Out[105]: 0.00025438479046726637
                In [106]: hypergeom.sf(3-1,100,5,10)
                # See Dog? You wnat num_common - 1
                Out[106]: 0.0066379128971171221
                # can we go back to drinking coffee now?
            """
            pval = hypergeom.sf(num_common - 1, num_universe, num_in_term, num_sampled)
            if pval <= pval_cutoff and num_common >= min_overlap:
                if label != None:
                    term.attrs["label"] = label
                term.attrs["hyper"] = OrderedDict(
                    [
                        ("pval", pval),
                        ("term_tested", len(terms)),
                        ("num_common", num_common),
                        ("num_universe", num_universe),
                        ("term_size", num_in_term),
                        ("num_terms", len(self)),
                        ("num_sampled", num_sampled),
                    ]
                )
                if bonferroni_correction == True:
                    # Right now this isn't true bonferroni, its only correcting for
                    # the number of terms that had term genes in it
                    if pval > pval_cutoff / len(terms):
                        term.attrs["hyper"]["bonferroni"] = False
                    else:
                        term.attrs["hyper"]["bonferroni"] = True
                term.attrs["pval"] = pval
                if include_genes == True:
                    term.attrs["hyper"]["genes"] = ",".join(
                        [x.id for x in term_genes.intersection(loci)]
                    )
                significant_terms.append(term)
        log.info(
            "\t\tFound {} was significant for {} terms", label, len(significant_terms)
        )
        if return_table == True:
            tbl = []
            for x in significant_terms:
                val = OrderedDict([("name", x.name), ("id", x.id)])
                val.update(x.attrs["hyper"])
                val.update(x.attrs)
                del val["hyper"]
                tbl.append(val)
            tbl = DataFrame.from_records(tbl)
            if label != None:
                tbl["label"] = label
            if len(tbl) > 0:
                tbl = tbl.sort_values(by="pval")
            return tbl
        else:
            return sorted(significant_terms, key=lambda x: x.attrs["pval"])

    # -----------------------------------------
    #       Internal Methods
    # -----------------------------------------

    def _initialize_tables(self):
        cur = self.m80.db.cursor()
        cur.execute(
            """
            CREATE TABLE IF NOT EXISTS terms (
                TID INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT UNIQUE,
                desc TEXT
            )"""
        )
        cur.execute(
            """
            CREATE TABLE IF NOT EXISTS term_loci (
                TID INTEGER, 
                LID INTEGER
            );"""
        )
        cur.execute(
            """
            CREATE TABLE IF NOT EXISTS term_attrs (
                TID INTEGER,
                key TEXT,
                val TEXT
            );
        """
        )
        cur.execute("CREATE INDEX IF NOT EXISTS term_index ON terms (name)")
        cur.execute("CREATE INDEX IF NOT EXISTS loci_index ON term_loci (TID,LID)")
        cur.execute("CREATE INDEX IF NOT EXISTS loci_attrs ON term_attrs (TID)")
        cur.execute("CREATE INDEX IF NOT EXISTS loci_attr_key ON term_attrs (key)")

    def _nuke_tables(self):
        cur = self.m80.db.cursor()
        cur.execute("DELETE FROM terms; DELETE FROM term_loci;")

    # -----------------------------------------
    #       Static Methods
    # -----------------------------------------

    @staticmethod
    def _parse_obo(obo_file) -> List[Term]:
        """
        Parses a Gene Ontology obo file
        Paramters
        ---------
        obo_file : filename
            The path the the obo file. You can download
            it here: http://geneontology.org/page/download-ontology
        
        Returns
        -------
        A list Terms with populated with name, desc, and attrs. OBO files
        do not have contain loci mappings, so terms will have no assigned 
        terms.
        """
        # Importing the obo information
        log.info(f"Importing OBO: {obo_file}")
        terms = []

        section = None
        name = None
        desc = None
        attrs = defaultdict(list)

        with rawFile(obo_file) as INOBO:
            for line in INOBO:
                line = line.strip()
                # Check to see if we have everything to add a term
                if line == "":
                    # NOTE: this needs to be separate if statements to avoid elif'ing on empty lines
                    if section == "[Term]":
                        terms.append(Term(name,desc=desc,attrs=attrs))
                        section,name,desc = None,None,None
                        attrs = defaultdict(list)
                # e.g. section line: [Term]
                elif line.startswith('[') and line.endswith(']'):
                    section = line
                else:
                    key,val = line.split(": ", maxsplit=1)
                    if " ! " in val:
                        val = val.split(" ! ", maxsplit=1)[0]
                    # A GO Term id is a name
                    if key == "id":
                        name = val
                    # A GO Term name is a desc
                    elif key == "name":
                        desc = val
                    else:
                        attrs[key].append(val)
        # Extract the children information
        log.info("Adding children attrs")
        rels = defaultdict(set)
        for t in terms:
            for is_a in t.attrs.get('is_a',[]):
                rels[t.name].add(is_a) 
        # Add children attrs to terms:
        for t in terms:
            if t.name in rels:
                t['children'] = list(rels[t.name])
        return terms

    @staticmethod
    def _parse_locus_term_map(
        locus_map_file,
        headers=True,
        go_col=1,
        id_col=0,
        sep="\t",
    ):
        # Importing locus map information, and cross referencing with obo information
        log.info(f"Importing Gene Map: {locus_map_file}")
        mapping = dict()
        locus = None
        term = None
        with rawFile(locus_map_file) as INMAP:
            if headers:
                INMAP.readline()
            for line in INMAP.readlines():
                if line.startswith("#") or line.startswith("!"):
                    continue
                row = line.strip("\n").split(sep)
                locus = row[id_col].split("_")[0].strip()
                term = row[go_col]
                # Make a map between loci and associated GO terms
                if term not in mapping:
                    mapping[term] = set([locus])
                else:
                    mapping[term].add(locus)
        return mapping

    # --------------------------------------------------
    #       factory methods
    # --------------------------------------------------

    @classmethod
    def from_terms(
        cls,
        name,
        terms,
        /,
        force=False,
        rootdir: Optional[str] = None,
    ):
        """
        Convenience function to create a Ontology from an iterable
        terms object.

        Parameters
        ----------
        terms : iterable of GOTerm objects
            Items to add to the ontology. The key being the name
            of the term and the items being the loci.
        name : str
            The name of the camoco object to be stored in the database.
        description : str
            A short message describing the dataset.
        loci : lp.Loci
            A RefGen object describing the genes in the dataset
        rootdir : str
            The base directory to store the files related to the dataset
            If not specified, 
        """
        if force:
            m80.delete("Ontology", name, rootdir=rootdir)
        # Do some checks
        if m80.exists("Ontology", name, rootdir=rootdir):
            raise ValueError(
                f"Ontology.{name} exists. Cannot use factory "
                f"methods on existing datasets."
            )
        ont = cls(name, rootdir=rootdir)
        try:
            # Store the loci information
            with ont.m80.db.bulk_transaction() as cur, ont.loci.m80.db.bulk_transaction() as lcur:
                for t in terms:
                    ont.add_term(t, cursor=cur, loci_cursor=lcur)
        except Exception as e:
            m80.delete("Ontology", name, rootdir=rootdir)
            raise e
        else:
            return ont

    @classmethod
    def from_obo(
        cls,
        name,
        obo_file,
        locus_map_file,
        loci,
        /,
        go_col=1,
        id_col=0,
        headers=True,
        force=False,
        rootdir: Optional[str] = None,
    ):
        """ 
        Convenience function for importing GOnt from obo files 
        Parameters
        ----------
        name : str
            The name of the camoco object to be stored in the database.
        obo_file : str
            Path to the obo file
        locus_map_file : str
            Path to the file which specifies what GO term each 
            locus is a part of.
        loci : lp.Loci
            A Loci instance to act as a reference for the Loci in the dataset
        go_col : int (default: 1)
            The index column for GO term in the locus_map_file
        id_col : int (default: 0)
            The index column for locus id in the locus_map_file
        headers : bool (default: True)
            A flag indicating whether or not there is a header line
            in the locus_map_file
        rootdir : str
            The base directory to store the files related to the dataset
            If not specified, 
        """
        if force:
            m80.delete("Ontology", name)
        if m80.exists("Ontology", name, rootdir=rootdir):
            raise ValueError(
                f"Ontology.{name} exists. Cannot use factory "
                f"methods on existing datasets."
            )
        ont = cls(name, rootdir=rootdir)
        # Define a helper function for extracting parents below
        def parents(terms, term):
            if "is_a" not in term.attrs:
                return None
            for t in term.attrs["is_a"]:
                yield terms[t]
                yield from parents(terms, terms[t]) 
        # Assemble ontology in try block, if fail, delete broken dataset
        try:
            terms = {term.name:term for term in ont._parse_obo(obo_file)}
            idmap = ont._parse_locus_term_map(locus_map_file)
            locimap = {l:loci[l] for l in set(chain(*idmap.values())) if l in loci} 
            missing_terms = set()
            # Start putting terms and loci together
            log.info("Populating Terms with Loci")
            for term_name,loci_names in idmap.items():
                if term_name not in terms:
                    missing_terms.add(term_name)
                    continue
                # Add the loci to the term
                term_loci = [locimap[l] for l in loci_names if l in loci]
                terms[term_name].loci.update(term_loci)
                # For each GO idmap, propagate loci to parent terms
                for parent in set(parents(terms,terms[term_name])):
                    parent.loci.update(term_loci)
            # Store the loci information
            with ont.m80.db.bulk_transaction() as cur, ont.loci.m80.db.bulk_transaction() as lcur:
                for t in terms.values():
                    ont.add_term(t, cursor=cur, loci_cursor=lcur)
        except Exception as e:
            m80.delete("Ontology", name, rootdir=rootdir)
            raise e
        else:
            return ont

