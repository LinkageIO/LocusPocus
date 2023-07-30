#!/usr/bin/python3


import math
import hashlib


import numpy as np

from typing import Any, Optional

from ..exceptions import StrandError, ChromosomeError, MissingLocusError
from .subloci import SubLoci
from .attrs import LocusAttrs


__all__ = ["Locus"]


class Locus:
    def __init__(
        self,
        chromosome: str,
        start: int,
        end: int,
        source: str = "locuspocus",
        feature_type: str = "locus",
        strand: str = "+",
        frame: int = None,
        name: str = None,
        # Extra locus stuff
        attrs: LocusAttrs = None,
        subloci: SubLoci = None,
    ):
        self.chromosome = str(chromosome)
        self.start = int(start)
        self.end = int(end)
        self.source = str(source)
        self.feature_type = str(feature_type)
        self.strand = str(strand)
        self.frame = frame
        self.name = name

        self.attrs = LocusAttrs(attrs)
        self.subloci = SubLoci(subloci)

    def __eq__(self, other):
        return hash(self) == hash(other) and self.attrs == other.attrs

    def __hash__(self):
        """
        Convert the locus to a hash, uses md5. The hash
        is computed using the *core* properties of the
        Locus, i.e. changing any attrs will not change
        the hash value.

        Parameters
        ----------
        None
        Returns
        -------
        int : md5 hash of locus
        """
        field_list = [
            str(x)
            for x in (
                self.chromosome,
                self.start,
                self.end,
                self.feature_type,
                self.strand,
                self.frame,
            )
        ]
        subloci_list = [str(hash(x)) for x in self.subloci]
        # Create a full string
        loc_string = "_".join(field_list + subloci_list)  # attr_list + subloci_list)
        digest = hashlib.md5(str.encode(loc_string)).hexdigest()
        return int(digest, base=16)

    def __len__(self):
        return abs(self.end - self.start) + 1

    def __lt__(self, locus):
        if self.chromosome == locus.chromosome:
            return self.start < locus.start
        else:
            return self.chromosome < locus.chromosome

    def __le__(self, locus):
        if (self.chromosome, self.coor) == (locus.chromosome, locus.coor):
            return True
        else:
            return self < locus

    def __ge__(self, locus):
        if (self.chromosome, self.coor) == (locus.chromosome, locus.coor):
            return True
        else:
            return self > locus

    def __gt__(self, locus):
        if self.chromosome == locus.chromosome:
            return self.start > locus.start
        else:
            return self.chromosome > locus.chromosome

    def __repr__(self):
        return (
            f"Locus("
            f"{self.chromosome},{self.start},{self.end},source={self.source},"
            f"feature_type='{self.feature_type}',"
            f"strand='{self.strand}',"
            f"frame='{self.frame}',"
            f"name='{self.name}',"
            f"attrs={self.attrs},"
            f"subloci=[{len(self.subloci)} subloci]"
            f")"
        )

    @property
    def stranded_start(self):
        if self.strand == "+":
            return min(self.coor)
        elif self.strand == "-":
            return max(self.coor)
        else:
            raise StrandError

    @property
    def stranded_end(self):
        if self.strand == "+":
            return max(self.coor)
        elif self.strand == "-":
            return min(self.coor)
        else:
            raise StrandError

    def __getitem__(self, item):
        return self.attrs[item]

    def __setitem__(self, key, val):
        self.attrs[key] = val

    def add_sublocus(
        self, locus: "Locus", find_parent=False, parent_attr: Optional[str] = "Parent"
    ):
        """
        Adds a sublocus to the current Locus. The added locus will be added
        according to its `parent_attr` keyword.
        """
        if not find_parent:
            self.subloci.add(locus)
        else:
            try:
                if locus[parent_attr] == self.name:
                    self.subloci.add(locus)
                else:
                    # Find the parent of the sublocus
                    parent = self.subloci.find(locus[parent_attr])
                    if parent is None:
                        raise MissingLocusError
                    parent.subloci.add(locus)
            except KeyError:
                raise KeyError(
                    f"Unable to resolve the key:{parent_attr} to find parent Locus"
                )

    def as_record(self):
        return (
            (
                self.chromosome,
                self.start,
                self.end,
                self.source,
                self.feature_type,
                self.strand,
                self.frame,
                self.name,
                hash(self),
            ),
            self.attrs,
        )

    def default_getitem(self, key, default=None) -> Any:
        """
        Returns the attr value of the Locus based on the key.
        similar to: Locus['Name'] -> "Locus1".
        If the attr key (e.g. 'Name') does not exist,
        a default can be specified.

        Parameters
        ----------
        key : str
            The key value of the attr
        default : Any
            If the key is not in the Locus attributes, return
            this value by default.
        """
        try:
            val = self.attrs[key]
        except KeyError:
            val = default
        finally:
            return val

    @property
    def coor(self):
        """
        Returns a tuple containing the start and end
        positions of the locus
        """
        return (self.start, self.end)

    def upstream(self, distance: int) -> int:
        """
        Calculates a base pair position 5' of the
        locus.

        Parameters
        ----------
        distance : int
            The distance upstream of the locus
        """
        if self.strand == "+":
            return max(0, self.start - distance)
        elif self.strand == "-":
            return self.end + distance

    def downstream(self, distance: int) -> int:
        """
        Calculates a base pair position 3' of the
        locus

        Parameters
        ----------
        distance : int
            The distance downstream of the locus
        """
        if self.strand == "+":
            return self.end + distance
        elif self.strand == "-":
            return max(0, self.start) - distance
        # return self.end + distance

    @property
    def center(self):
        """
        Calculates the center base pair position of
        the locus.

        NOTE: If the locus has an odd length, a 'half-bp'
              will be returned.
              E.g: Locus('1',100,200).center ==

        Returns
        -------
        The center position
        """
        return self.start + len(self) / 2

    def distance(self, locus):
        """
        Return the number of base pairs between two loci.
        NOTE: this excludes the start/end bases of the loci.

        Locus A             Locus B
        ==========---------=========
        1       10         20      30

        There are 9 bases between 10 and 20 (excluding positions
        10 and 20 themselves because they are included in the
        loci).

        If the loci are on different chromosomes, return np.inf

        Parameters
        ----------
        locus: Locus Object
            A second locus object to calculate distance.

        Returns
        -------
        int: the number of bp between the loci
        np.inf: if on different chromosomes

        """
        if self.chromosome != locus.chromosome:
            distance = np.inf
        else:
            x, y = sorted([self, locus])
            distance = y.start - x.end - 1
        return distance

    def center_distance(self, locus):
        """
        Return the distance between the center of two loci.
        If the loci are on different chromosomes, return np.inf.

        Parameters
        ----------
        locus : Locus Object
            A second locus object used to calculate a distance.

        Returns
        -------
        int : the distance between the center of two loci.
        np.inf: if on different chromosomes

        """
        if self.chromosome != locus.chromosome:
            distance = np.inf
        else:
            distance = math.floor(abs(self.center - locus.center))
        return distance

    def combine(self, locus):
        """
        Returns a new Locus with start and stop boundaries
        that contain both of the input loci. Both input loci are
        added to the subloci of the new Locus.

        NOTE: this ignores strand, the resultant Locus is just
              a container for the input loci.

        ___________Ascii Example__________________________

              Locus A          Locus B
        ------=========-------=============---------------

        A.combine(B)
        ------=============================---------------
              subloci=[A,B]
        """
        if self.chromosome != locus.chromosome:
            raise ChromosomeError("Input Chromosomes do not match")
        x, y = sorted([self, locus])
        start = x.start
        end = y.end
        return Locus(self.chromosome, start, end, subloci=[self, locus])

    def as_tree(self, parent=None):  # pragma: no cover
        from anytree import Node, RenderTree

        root = Node(f"{self.feature_type}:{self.name}", parent=parent)
        for c in self.subloci:
            node = c.as_tree(parent=root)
        if parent is None:
            for pre, _, node in RenderTree(root):
                print("%s%s" % (pre, node.name))
        return root

    def __str__(self):
        return repr(self)

    # --------------------------------
    #   Factory Methods
    # --------------------------------
    @classmethod
    def from_gff_line(
        cls,
        line,
        /,
        ID_attr: str = "ID",
        parent_attr: str = "Parent",
        attr_split: str = "=",
    ) -> "Locus":
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
        ) = line.strip().split("\t", maxsplit=8)
        # Cast data into appropriate types
        strand = None if strand == "." else strand
        frame = None if frame == "." else int(frame)
        # Get the attributes
        attributes = dict(
            [
                (field.strip().split(attr_split))
                for field in attributes.strip(";").split(";")
            ]
        )
        # Store the score in the attrs if it exists
        if score != ".":
            attributes["score"] = float(score)
        # Parse out the Identifier
        if ID_attr in attributes:
            name = attributes[ID_attr]
        else:
            name = None
        l = cls(
            chromosome,
            start,
            end,
            source=source,
            feature_type=feature,
            strand=strand,
            frame=frame,
            name=name,
            attrs=attributes,
        )
        return l
