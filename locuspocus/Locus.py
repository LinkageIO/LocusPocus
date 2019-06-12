#!/usr/bin/python3
from collections import defaultdict
from dataclasses import dataclass,field,InitVar,FrozenInstanceError
from itertools import chain
from typing import Union, Any, List, Optional, cast, Callable, Iterable, Dict

from .Exceptions import StrandError, ChromosomeError

import re
import math
import locuspocus
import dataclasses
import uuid

import pandas as pd
import numpy as np

from locuspocus.Loci import Loci


__all__ = ['Locus']

class Locus:
    
    _ref = Loci()
    
    def __init__(
        self,
        chromosome: str,
        start: int,
        end: int,

        source: str = 'locuspocus',
        feature_type: str = 'locus',
        strand: str = '+',
        frame: int = None,
        name: str = None,

        # Extra locus stuff
        attrs: dict = None,
        parent = None,
        children = None
    ):
        if (chromosome,start,end) == (None,None,None):
            # This allows Locus objects to be created by class methods 
            self._LID = name
        else:
            # Add the info to the database
            LID = self._ref.add_locus(
                chromosome,
                start,
                end,
                source,
                feature_type,
                strand,
                frame,
                name,
                attrs,
                parent,
                children
            )
            self._LID = LID

    def __del__(self):
        '''
            Deletes a locus object, inlcuding any traces 
            of it from the Loci db
        '''
        self._ref.del_locus(self)

    def _core_property(self,key):
        val, = self._ref._db.cursor().execute(f'''
            SELECT {key} FROM loci WHERE LID = ?
        ''',(self._LID,)
        ).fetchone()
        return val

    @property
    def chromosome(self):
        return self._core_property('chromosome')
        
    @property
    def start(self):
        return self._core_property('start')

    @property
    def end(self):
        return self._core_property('end')

    @property
    def feature_type(self):
        return self._core_property('feature_type')

    @property
    def strand(self):
        return self._core_property('strand')

    @property
    def frame(self):
        return self._core_property('frame')

    @property
    def name(self):
        return self._core_property('name') 

    @property
    def parent(self):
        parent_lid, = self._ref._db.cursor('''
            SELECT parent FROM relationships WHERE child = ?
        ''',(self._LID))
        return self.from_LID(parent_LID)

    @parent.setter
    def parent(self,val):
        raise NotImplementedError

    @property
    def children(self):
        raise NotImplementedError

    @children.setter
    def children(self,val):
        raise NotImplementedError


    def __len__(self):
        return abs(self.end - self.start) + 1

    def __eq__(self,locus):
        try:
            assert self.chromosome == locus.chromosome 
            assert self.start == locus.start 
            assert self.end == locus.end 
            assert self.source == locus.source 
            assert self.feature_type == locus.feature_type 
            assert self.strand == locus.strand 
            assert self.frame == locus.frame 
            assert self.name == locus.name 
            for key,val in self.attrs.items():
                assert val == locus.attrs[key]
            for x,y in zip(sorted(self.subloci),sorted(locus.subloci)):
                assert x == y
            return True
        except (AssertionError,KeyError):
            return False

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

    @classmethod
    def from_LID(cls,LID,ref=None):
        # Create a locus using sentinels
        l = cls(None,None,None,name=LID)
        # assign ref if needed
        if ref is not None:
            l._ref = ref
        return l

