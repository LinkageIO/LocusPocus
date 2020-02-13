#!/usr/bin/python3
from collections import defaultdict
from dataclasses import dataclass,field,InitVar,FrozenInstanceError
from itertools import chain
from typing import (
    Union, Any, List, Optional, 
    cast, Callable, Iterable
)

from .exceptions import StrandError, ChromosomeError

import re
import math
import locuspocus
import dataclasses

import pandas as pd
import numpy as np

__all__ = ['Locus']

class SubLoci():
    # A restricted list interface to subloci
    def __init__(self,loci=None):
        self._loci = loci

    @property
    def empty(self):
        if self._loci is None:
            return True
        else:
            return False

    def __len__(self):
        if self.empty:
            return 0
        else:
            return len(self._loci)

    def __eq__(self, other):
        if len(self) != len(other):
            return False
        else:
            return sorted(self) == sorted(other)

    def __iter__(self):
        if self.empty:
            return (x for x in [])
        return (x for x in self._loci)

    def add(self,locus):
        if self.empty:
            self._loci = []
        self._loci.append(locus)

    def __getitem__(self,index):
        if self.empty:
            raise IndexError
        return self._loci[index]

    def __len__(self):
        if self.empty:
            return 0
        return len(self._loci)

    def __repr__(self):
        if len(self) == 0 or self.empty:
            return "[]"
        return f"[{len(self)} subloci]"
    

class LocusAttrs():
    # a restricted dict interface to attributes
    def __init__(self,attrs=None):
        self._attrs = attrs

    def __len__(self):
        if self.empty:
            return 0
        else:
            return len(self._attrs)

    def __eq__(self, other):
        if self.empty and other.empty:
            return True
        elif len(self) != len(other):
            # Short circuit on length
            return False
        else:
            return sorted(self.items()) == sorted(other.items())

    @property
    def empty(self):
        if self._attrs is None:
            return True
        else:
            return False

    def keys(self):
        if self.empty:
            return []
        else:
            return self._attrs.keys()

    def values(self):
        if self.empty:
            return []
        else:
            return self._attrs.values()

    def items(self):
        if self.empty:
            return {}
        else:
            return self._attrs.items()

    def __contains__(self,key):
        if self.empty:
            return False
        return key in self._attrs

    def __getitem__(self,key):
        if self.empty:
            raise KeyError()
        return self._attrs[key]

    def __setitem__(self,key,val):
        if self.empty:
            self._attrs = {}
        self._attrs[key] = val

    def __repr__(self):
        if self.empty:
            return repr({})
        return repr(self._attrs)

@dataclass()
class Locus:
    chromosome: str
    start: int 
    end: int

    source: str = 'locuspocus'
    feature_type: str = 'locus'
    strand: str = '+'
    frame: int = None
    name: str = None

    # Extra locus stuff
    attrs: LocusAttrs = None 
    subloci: SubLoci  = None 

    def __post_init__(self):
        self.attrs = LocusAttrs(self.attrs)
        self.subloci = SubLoci(self.subloci)

    def __eq__(self,other):
        if (self.chromosome == other.chromosome
            and self.start == other.start
            and self.end == other.end
            and self.source == other.source
            and self.feature_type == other.feature_type
            and self.strand == other.strand
            and self.frame == other.frame
            and self.name == other.name
            and self.attrs == other.attrs
            and self.subloci == other.subloci
        ):
            return True
        else:
            return False


    def __len__(self):
        return abs(self.end - self.start) + 1

    def __lt__(self,locus):
        if self.chromosome == locus.chromosome:
            return self.start < locus.start
        else:
            return self.chromosome < locus.chromosome

    def __le__(self,locus):
        if (self.chromosome,self.coor) == (locus.chromosome,locus.coor):
            return True
        else:
            return self < locus

    def __ge__(self,locus):
        if (self.chromosome,self.coor) == (locus.chromosome,locus.coor):
            return True
        else:
            return self > locus

    def __gt__(self,locus):
        if self.chromosome == locus.chromosome:
            return self.start > locus.start
        else:
            return self.chromosome > locus.chromosome

    @property
    def stranded_start(self):
        if self.strand == '+':
            return min(self.coor)
        elif self.strand == '-':
            return max(self.coor)
        else:
            raise StrandError

    @property
    def stranded_end(self):
        if self.strand == '+':
            return max(self.coor)
        elif self.strand == '-':
            return min(self.coor)
        else:
            raise StrandError

    def __getitem__(self,item):
        return self.attrs[item]

    def __setitem__(self,key,val):
        self.attrs[key] = val

    def add_sublocus(self,locus):
        self.subloci.add(locus)

    def as_record(self):
        return ((
            self.chromosome,
            self.start,
            self.end,
            self.source,
            self.feature_type,
            self.strand,
            self.frame,
            self.name,
        ),
            self.attrs
        )

    def default_getitem(self,key,default=None) -> Any:
        '''
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
        '''
        try:
            val = self.attrs[key]
        except KeyError as e:
            val = default
        finally:
            return val

    @property
    def coor(self):
        '''
        Returns a tuple containing the start and end
        positions of the locus
        '''
        return (self.start,self.end)

    def upstream(self,distance: int) -> int:
        '''
        Calculates a base pair position 5' of the 
        locus.

        Parameters
        ----------
        distance : int
            The distance upstream of the locus
        '''
        if self.strand == '+':
            return max(0,self.start - distance)
        elif self.strand == '-':
            return self.end + distance

    def downstream(self,distance: int) -> int:
        '''
        Calculates a base pair position 3' of the 
        locus
        
        Parameters
        ----------
        distance : int
            The distance downstream of the locus
        '''
        if self.strand == '+':
            return self.end + distance
        elif self.strand == '-':
            return max(0,self.start) - distance
        #return self.end + distance

    @property
    def center(self):
        '''
        Calculates the center base pair position of
        the locus.

        NOTE: If the locus has an odd length, a 'half-bp'
              will be returned. 
              E.g: Locus('1',100,200).center == 

        Returns
        -------
        The center position
        '''
        return self.start + len(self)/2

    def distance(self,locus):
        '''
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

        '''
        if self.chromosome != locus.chromosome:
            distance = np.inf
        else:
            x,y = sorted([self,locus])
            distance = y.start - x.end - 1
        return distance

    def center_distance(self, locus):
        '''
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

        '''
        if self.chromosome != locus.chromosome:
            distance = np.inf
        else:
            distance = math.floor(abs(self.center - locus.center))
        return distance

    def combine(self, locus):
        '''
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
        '''
        if self.chromosome != locus.chromosome:
            raise ChromosomeError('Input Chromosomes do not match') 
        x,y = sorted([self,locus])
        start = x.start
        end = y.end
        return Locus(self.chromosome,start,end,subloci=[self,locus])

    def as_tree(self,parent=None): #pragma: no cover
        from anytree import Node, RenderTree
        root = Node(f'{self.feature_type}:{self.name}',parent=parent)
        for c in self.subloci:
            node = c.as_tree(parent=root)
        if parent is None:
            for pre, _, node in RenderTree(root):
                print("%s%s" % (pre, node.name))
        return root

    def __str__(self):
        return repr(self)
