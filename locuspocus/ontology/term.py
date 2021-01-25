#!/usr/bin/python3

import numpy as np


class Term:
    """
    A Term is a named set of loci.

    NOTE: this is different that a Loci obect which is a
          named set of **reference** loci. Loci within a term
          are somehow related ouside the context of the whole
          genome, for instance, in some biological function.

    Parameters
    ----------
    name : unique identifier
    desc: short description
    loci : iterable of loci objects that are related
    attrs : dictionary of term attributes

    Returns
    -------
    A Term Object

    """

    def __init__(self, name, desc="", loci=None, attrs=None):
        self.name = name
        self.desc = desc
        # Manage the loci
        self.loci = set()
        if loci is not None:
            [self.add_locus(l) for l in loci]
        # Manage the attrs 
        self.attrs = dict() 
        if attrs is not None:
            for k,v in attrs.items():
                if isinstance(v,(int,float,str)):
                    self.attrs[k] = [v]
                elif isinstance(v,list) and all([isinstance(_,(int,float,str)) for _ in v]):
                    self.attrs[k] = v
                else:
                    raise TypeError("Attrs vals must be one of OR a list of: int,float,str")

    def __len__(self):
        """
        Returns the number of loci in the term.
        """
        return len(self.loci)

    def __getitem__(self, key):
        val = self.attrs[key]
        if len(val) == 1:
            val = val[0]
        return val

    def __setitem__(self, key, val):
        self.attrs[key] = [val]

    def attr_append(self, key, val):
        if key not in self.attrs:
            self.attrs[key] = []
        self.attrs[key].append(val)

    def add_locus(self, locus):
        """
        Adds a locus to the Term.
        """
        self.loci.add(locus)

    def nearby_loci(self, locus, max_distance):
        """
        Returns Term loci within max_distance
        base pairs away from the input locus.

        Parameters
        ----------
        locus: Locus Object
            The input locus from which to identify
            nearby Term loci.
        max_distance: int
            The maximum distance considered to be
            'nearby'

        Returns
        -------
        a generator that yields nearby loci

        """
        for l in self.loci:
            if locus.distance(l) <= max_distance:
                yield l

    def effective_loci(self, max_distance):
        """
        Collapse down loci within max_distance into
        'effective' loci. Looks like:

        ____________Ascii Example____________________
               Locus1         Locus2
        -------========-------=========--------------
              50     100     150    200
        max_distance=100

               Locus3
        -------========================--------------
            subloci:
               ========       =========


        Parameters
        ----------
        window_size : int
            The maximum distance two loci need to be to
            not be collapsed into an effective locus
        """
        loci = sorted(self.loci)
        collapsed = [loci.pop(0)]
        for locus in loci:
            tail = collapsed[-1]
            # if they have overlapping windows, collapse
            if tail.distance(locus) <= max_distance:
                collapsed[-1] = tail.combine(locus)
            else:
                collapsed.append(locus)
        print(
            f"Term({self.name}): {len(self.loci)} Loci -> "
            f"{len(collapsed)} effective Loci "
            f"within {max_distance} bp"
        )
        return collapsed

    def strongest_loci(
        self,
        attr,
        window_size=None,
        lowest=True
    ):
        """
        Collapses down loci that have overlapping windows,
        then returns the locus with the strongest 'attr'
        value. Looks like:

        Locus1:    |--------o-------| (attr: 7)
        Locus2:        |--------o--------| (attr: 2)
        Locus3:                             |--------o--------| (attr: 8)
        Strongest: |--------o-------|       |--------o--------|

        Legend: '|' : Window edge, used to collapse
                'o' : 'Locus' edge (SNPs in this case)

        Parameters
        ----------
        attr : str
            The locus attribute to use to determine the 'strongest'
        window_size : int (default: None)
            If not None, maps a new window size to each locus.
        lowest: bool (default: True)
            When sorting by attr, lowest is strongest (i.e. p-vals)
        """
        raise NotImplementedError
        is_reverse = not lowest
        return [
            # sort by attr and take first item
            sorted(
                locus.sub_loci,
                key=lambda x: float(x.default_getitem(attr, np.inf)),
                reverse=is_reverse,
            )[0]
            for locus in self.effective_loci(window_size=window_size)
        ]
