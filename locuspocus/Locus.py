#!/usr/bin/python3
from collections import defaultdict
from dataclasses import dataclass,field,InitVar,FrozenInstanceError
from itertools import chain
from typing import Union, Any, List, Optional, cast, Callable, Iterable

from .Exceptions import StrandError

import re
import math
import hashlib
import locuspocus
import dataclasses

import pandas as pd
import numpy as np

class SubLoci(object):
    '''
        This class can either be an actual list of subloci OR
        a list of LIDS and a reference to the database where the 
        subloci can be accessed as they are needed.
    '''
    def __init__(self, subloci, refloci=None):
        self._subloci = subloci if subloci is not None else []
        if refloci is None:
            # If there is not db ref, make sure all items are Locuses
            assert all([isinstance(x,Locus) for x in self._subloci])
            self.db_backed = False
        else:
            assert all([isinstance(x,int) for x in self._subloci])
            self._refloci = refloci
            self.db_backed = True

    def __getitem__(self, key):
        if not self.db_backed :
            return self._subloci[key]
        else:
            return self._refloci._get_locus_by_LID(self._subloci[key])

    def __iter__(self):
        if not self.db_backed:
            yield from (x for x in self._subloci)
        else:
            yield from (self._refloci._get_locus_by_LID(i) for i in self._subloci)

    def __repr__(self): #pragma: no cover
        return repr(list(self))

    def __len__(self):
        return len(self._subloci)

    def append(self,x):
        if self.db_backed:
            raise ValueError("cannot change subloci of a db backed Locus")
        else:
            assert isinstance(x,Locus)
            self._subloci.append(x)

class LociAttrs(object):
    '''
        This class can either be a normal dictionary with attributes OR
        a reference to a RefLoci database and a LID for the locus.
    '''
    def __init__(self, attrs, refloci=None, LID=None):
        if refloci is None and LID is None:
            if attrs is None:
                self._attrs = {}
            else:
                assert isinstance(attrs,dict)
                self._attrs = attrs
            self.db_backed = False
        elif refloci is not None and LID is not None:
            self._LID = LID
            self._refloci = refloci
            self.db_backed = True

    def _db_getitem(self,key):
        if not self.db_backed:
            raise ValueError('The locus is not backed by a database')
        cur = self._refloci._db.cursor()
        results = cur.execute(
            '''SELECT val FROM loci_attrs WHERE LID = ? AND key = ? ''',
            (self._LID,key)
        ).fetchone()
        if results is None:
            raise KeyError(f'{key} not in attrs')
        else:
            return results[0]

    def __getitem__(self,key):
        if not self.db_backed:
            return self._attrs[key]
        else:
            return self._db_getitem(key)

    def __setitem__(self,key,val):
        raise ValueError('Loci attrs cannot be changed')

    def keys(self):
        if not self.db_backed:
            return self._attrs.keys()
        else:
            cur = self._refloci._db.cursor()
            keys = cur.execute(
                '''SELECT key FROM loci_attrs WHERE LID = ?''',
                (self._LID,)
            ).fetchall()
            return (k[0] for k in keys)

    def values(self): 
        if not self.db_backed:
            return self._attrs.values()
        else:
            cur = self._refloci._db.cursor()
            vals = cur.execute(
                '''SELECT val FROM loci_attrs WHERE LID = ?''',
                (self._LID,)
            ).fetchall()
            return (v[0] for v in vals)

    def items(self):
        if not self.db_backed:
            return self._attrs.items()
        else:
            return (
                (k,v) for k,v in zip(self.keys(),self.values())        
            )

    def __contains__(self,key):
        if not self.db_backed:
             return key in self._attrs
        else:
            try:
                val = self._db_getitem(key)
                return True
            except KeyError as e:
                return False

    def __repr__(self):
        return str(dict(self.items()))

@dataclass()
class Locus:
    chromosome: str
    start: int 
    end: int

    source: str = 'locuspocus'
    feature_type: str = None
    strand: str = '+'
    frame: int = None

    # Extra locus stuff
    name: str = None
    attrs: InitVar[attrs] = field(default=None)
    subloci: InitVar[subloci] = field(default=None) 
    refloci: InitVar[str] = field(default=None)
    _frozen: bool = field(default=False,repr=False)
    _LID: int = field(default=None,repr=False)

    def __post_init__(self, attrs, subloci, refloci):
        self.refloci = refloci
        # Handle the "smart" things
        self.subloci = SubLoci(subloci,refloci)
        self.attrs = LociAttrs(attrs,refloci,self.LID)
        # Freeze the object
        self._frozen = True

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

    def __hash__(self):
        """
            Convert the locus to a hash, uses md5.

            Parameters
            ----------
            None

            Returns
            -------
            str : md5 hash of locus

            Notes
            -----
            The fields used for this computation are: 
            (chromosome, start, end, source, feature_type,
            strand, frame, name). Collisions will happen
            if two loci only differ by attrs or subloci.
            This is **NOT** a uuid.

        """
        loc_string  = "_".join([str(x) for x in (
            self.chromosome,
            self.start,
            self.end,
            self.source,
            self.feature_type,
            self.strand,
            self.frame,
            self.name
        )])
        digest = hashlib.md5(str.encode(loc_string)).hexdigest()
        return int(digest, base=16)

    @property
    def LID(self):
        return self._LID

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
        if self.attrs is not None: 
            return self.attrs[item]

    def __setitem__(self,key,val):
        raise FrozenInstanceError("Cannot change attrs of Locus")

    def __setattr__(self,key,val):
        if self._frozen is True:
            raise FrozenInstanceError("Cannot change attrs of Locus")
        else:
            super().__setattr__(key,val)

    def add_sublocus(self,locus):
        self.subloci.append(locus)

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
            hash(self)
        ),
            self.attrs
        )

    def as_dict(self) -> dict:
        '''
        Returns the Locus as a dictionary
        '''
        return dataclasses.asdict(self) 

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

    def center_distance(self, locus):
        '''
        Return the absolute distance between the center of two loci. 
        If the loci are on different chromosomes, return np.inf.

        Parameters
        ----------
        locus : Locus Object
            A second locus object used to calculate a distance.

        Returns
        -------
        int : the distance between the center of two loci.
        '''
        if self.chromosome != locus.chromosome:
            distance = np.inf
        else:
            distance = math.floor(abs(self.center - locus.center))
        return distance

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
