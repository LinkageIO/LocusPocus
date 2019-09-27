import re
import math
import locuspocus
import dataclasses
import uuid

import numpy as np

from typing import Union, Any, List, Optional, cast, Callable, Iterable, Dict
from locuspocus.exceptions import StrandError,ChromosomeError,LocusError

class LocusMixin(object):

    def _core_property(self,key):
        val, = self._db.cursor().execute(f'''
            SELECT {key} FROM loci WHERE UUID = ?
        ''',(self._UUID,)
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
    def source(self):
        return self._core_property('source')

    @property
    def score(self):
        try:
            return self['score']
        except KeyError:
            return None

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
    def attrs(self):
        class attrView:
            def __init__(self,UUID,db):
                self._UUID = UUID
                self._db = db

            def keys(self):
                return [ x[0] for x in 
                    self._db.cursor().execute('''
                        SELECT key FROM loci_attrs
                        WHERE UUID = ?;
                    ''',(self._UUID,)).fetchall()
                ]
            def values(self):
                return [ x[0] for x in 
                    self._db.cursor().execute('''
                        SELECT val FROM loci_attrs
                        WHERE UUID = ?;
                    ''',(self._UUID,)).fetchall()
                ]
            def items(self):
                return [ (x[0],x[1]) for x in 
                    self._db.cursor().execute('''
                        SELECT key,val FROM loci_attrs
                        WHERE UUID = ?;
                    ''',(self._UUID,)).fetchall()
                ]
            def __contains__(self,key):
                count, = self._db.cursor().execute('''
                    SELECT COUNT(*) FROM loci_attrs
                    WHERE UUID = ? AND key = ?;
                ''',(self._UUID,key)).fetchone()
                return True if count == 1 else False

        return attrView(self._UUID,self._db)

    def __getitem__(self,key):
        '''
        Retrieving attributes that are not the core properties requires
        a string key. 

        Example:
        >>> l = lp.Locus('1',1,100,attrs={'foo':'bar'})
        >>> l['foo']
        'bar'

        '''
        try:
            val,val_type = self._db.cursor().execute('''
                SELECT val,type FROM loci_attrs 
                WHERE UUID = ?
                AND key = ?
            ''',(self._UUID,key)).fetchone()
        except (ValueError,TypeError) as e:
           raise KeyError(f'{key} not in Locus attrs')
        if val_type == 'int':
            val = int(val)
        elif val_type == 'float':
            val = float(val)
        elif val_type == 'bool':
            val = bool(val)
        elif val_type == 'str':
            val = str(val)
        return val

    def __setitem__(self,key,val):
        if not isinstance(key,str):
            raise ValueError('Locus attribute key must be a str')
        if not isinstance(val,(str,int,float,bool)):
            raise ValueError(
                'Locus attribute value must '
                'be one of: (str,int,float,bool)'
            )
        val_type = guess_type(val)
        self._db.cursor().execute('''
            INSERT OR REPLACE INTO loci_attrs
            (UUID,key,val,type)
            VALUES
            (?,?,?,?)
        ''',(self._UUID,key,val,val_type))

    # Tree methods
    def __detach(self):
        ' detach self from a parent locus'
        cur = self._db.cursor()
        cur.execute('''
            DELETE FROM relationships
            WHERE child = ?;
        ''',(self._UUID,))

    @property
    def parent(self):
        try:
            parent_UUID, = self._db.cursor().execute('''
                SELECT parent FROM relationships WHERE child = ?
            ''',(self._UUID,)).fetchone()
            parent = self._from_UUID(parent_UUID)
        except TypeError:
            parent =  None
        return parent

    @parent.setter
    def parent(self,parent):
        if parent is None:
            self.__detach()
        else:
            try:
                with self._db:
                    cur = self._db.cursor()
                    # detach self from tree
                    self.__detach()
                    cur.execute('''
                        INSERT OR REPLACE INTO relationships
                        (parent,child) 
                        VALUES
                        (?,?)
                    ''',(parent._UUID,self._UUID))
            except AttributeError as e:
                raise LocusError('Parent must be type: Locus') 

    @property
    def children(self):
        children = [
            self._from_UUID(x[0]) for x in \
            self._db.cursor().execute('''
                SELECT child FROM relationships WHERE parent = ?
            ''',(self._UUID,)
            )
        ]
        return tuple(children)

    @children.setter
    def children(self,children):
        # detach current children
        for child in self.children:
            child.__detach()
        # attach new children
        if children is not None:
            for child in children:
                child.parent = self
        
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
            for k,v in self.attrs.items():
                assert locus[k] == v
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
            val = self[key]
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
        added to the children of the new Locus.

        NOTE: this ignores strand, the resultant Locus is just
              a container for the input loci.

        ___________Ascii Example__________________________
             
              Locus A          Locus B
        ------=========-------=============---------------
        
        A.combine(B)
        ------=============================---------------
              children=[A,B]
        '''
        if self.chromosome != locus.chromosome:
            raise ChromosomeError('Input Chromosomes do not match') 
        x,y = sorted([self,locus])
        start = x.start
        end = y.end
        return self.create(self.chromosome,start,end,children=[self,locus])

    def as_tree(self,parent=None): #pragma: no cover
        from anytree import Node, RenderTree
        root = Node(f'{self.feature_type}:{self.name}',parent=parent)
        for c in self.children:
            node = c.as_tree(parent=root)
        if parent is None:
            for pre, _, node in RenderTree(root):
                print("%s%s" % (pre, node.name))
        return root

    def __str__(self):
        return repr(self)

    def __repr__(self):
        if self.name is not None:
            string = f"{self.feature_type}('{self.name}')"
        else:
            string =  (
                f"{self.feature_type}({self.chromosome},{self.start},{self.end},"
                f"source='{self.source}',"
                f"strand='{self.strand}',"
                f"frame={self.frame},"
                f"name='{self.name}'"
                ')'
            )
        return string

def guess_type(object):
    '''
        Guess the type of object from the class attribute
    '''
    # retrieve a list of classes
    classes = re.match(
        "<class '(.+)'>",
        str(object.__class__)
    ).groups()[0].split('.')
    # Return the most specific one
    return classes[-1]
