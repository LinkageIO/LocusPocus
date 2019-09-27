#!/usr/bin/python3

import logging
import numpy as np

from dataclasses import dataclass,field,InitVar,FrozenInstanceError

@dataclass
class Term:
    '''
        A Term is a named set of loci.

        NOTE: this is different that a RefLoci obect which is a 
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

    '''

    name: str = field(default=None)
    desc: str = field(default=None)
    loci: set = field(default_factory=list)
    attrs: dict = field(default_factory=dict)
    parent: None = field(default=None)
    _frozen: bool = field(default=False,repr=False)
    _TID: int = field(default=None,repr=False)

    def __post_init__(self):
        # Force this coversion
        self.loci = set(self.loci)
        self._frozen = True

    def __len__(self):
        '''
        Returns the number of loci in the term.
        '''
        return len(self.loci)

    def __getitem__(self,key):
        return self.attrs[key]

    def __setitem__(self,key,val):
        if self._frozen is True and key != '_frozen':
            raise FrozenInstanceError("Cannot change attrs of Locus")
        else:
            super().__setattr__(key,val)

    def add_locus(self, locus):
        '''
        Adds a locus to the Term.
        '''
        self.loci.add(locus)

    def nearby_loci(self, locus, max_distance):
        '''
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

        '''
        for l in self.loci:
            if locus.distance(l) <= max_distance: 
                yield l
    
    def effective_loci(self, max_distance):
        '''
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
        '''
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
            f'Term({self.name}): {len(self.loci)} Loci -> '
            f'{len(collapsed)} effective Loci '
            f'within {max_distance} bp'
        )
        return collapsed

   #def strongest_loci(self, attr, window_size=None, lowest=True):
   #    '''
   #        Collapses down loci that have overlapping windows,
   #        then returns the locus with the strongest 'attr'
   #        value. Looks like:

   #        Locus1:    |--------o-------| (attr: 7)
   #        Locus2:        |--------o--------| (attr: 2)
   #        Locus3:                             |--------o--------| (attr: 8)
   #        Strongest: |--------o-------|       |--------o--------|

   #        Legend: '|' : Window edge, used to collapse
   #                'o' : 'Locus' edge (SNPs in this case)

   #        Parameters
   #        ----------
   #        attr : str
   #            The locus attribute to use to determine the 'strongest'
   #        window_size : int (default: None)
   #            If not None, maps a new window size to each locus.      
   #        lowest: bool (default: True)
   #            When sorting by attr, lowest is strongest (i.e. p-vals) 
   #    '''
   #    is_reverse = not lowest
   #    return [
   #        # sort by attr and take first item
   #        sorted(
   #            locus.sub_loci,
   #            key=lambda x: float(x.default_getitem(attr,np.inf)),
   #            reverse=is_reverse
   #        )[0] for locus in self.effective_loci(window_size=window_size)
   #    ]

    def copy(
        self, 
        name=None, 
        desc=None,
        loci=None,
        attrs=None
    ):
        '''
        Creates a copy of a term with the option to 
        expand loci and attrs. 

        Parameters
        ----------
        name : str
            A required name for the new term.
        desc : str
            An optional short description for the term.
        loci : iterable of co.Loci objects
            These loci will be added to the Term object
            in addition to the loci objects that were
            in the original Term.
        attrs : key value pairs
            Additional key value pairs will be added 
            as attributes to the term object.

        Returns
        -------
        A Term object.
        '''
        if name is None:
            name = self.name
        if loci == None:
            loci = self.loci.copy()
        if attrs is None:
            attrs = dict()
        loci = self.loci.extend(loci)
        new_attrs = self.attrs.copy()
        new_attrs.update(**attrs)
        copy = Term(
            name,
            desc=desc,
            loci=loci,
            attrs=new_attrs
        )
        return copy

