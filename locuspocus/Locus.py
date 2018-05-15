#!/usr/bin/python3
from collections import defaultdict
from itertools import chain

import hashlib
import re
import pandas as pd

class Locus(object):
    def __init__(self, chrom, start, end=None,
                 id=None, global_order=0, feature_order=0,
                 feature_type=None, window=0, sub_loci=None, **kwargs):
        '''
            A genetic coordinate class.

            Looks like:

              Window
              ___|____
             /        \
            |----------=========---------|
            ^          ^       ^         ^
             \          \       \         \
              Upstream   start   end       downstream

        '''
        # Intelligently assign an ID, which is not required
        if id is None or id.startswith('<None>') or id == '.':
            self._id = None
        else:
            self._id = id
        # Chromosomes are strings
        self.chrom = str(chrom)
        self._feature_type = str(feature_type) if feature_type is not None else str(None)
        # Postitions are integers
        self._start = int(start)
        self._end = int(end) if end is not None else int(start)
        self._global_order = int(global_order)
        self._feature_order = int(feature_order)

        # Implement an optional window around the start and stop
        # This is used to collapse SNPs and to perform the
        # upstream and downstream methods
        self.window = int(window)
        # Keep a dictionary for special locus attributes
        self.attr = kwargs

        # Loci can also have sub loci (for exons, snps, etc)
        self.sub_loci = set(sub_loci) if sub_loci is not None else set()
        if len(self.sub_loci) == 0:
            self.sub_loci.add(self)
        #  Warn us if something seems off
        if self._start > self._end:
            raise ValueError(
                "Wonky start and stop positions for: {}".format(self)
            )

    def as_dict(self):
        '''
            Return the Locus information as a dictionary

            Parameters
            ----------
            None!

            Returns
            -------
            a python dictionary

        '''
        a_dict = {
            'id'    : self.id,
            'chrom' : self.chrom,
            'start' : self.start,
            'end'   : self.end,
            'global_order': self.global_order,
            'feature_order': self.feature_order,
            'feature_type': self.feature_type
        }
        a_dict.update(self.attr)
        return a_dict

    @property
    def id(self):
        '''
            Return the Locus id

            Parameters
            ----------
            None

            Returns
            -------
            id : string
        '''
        if self._id is None:
            return '''<{}>{}:{}-{};[{},{},{}]'''.format(
                self._id, self.chrom,
                self.start, self.end,
                self.global_order, self.feature_order, self.feature_type
            )
        else:
            return self._id

    def update(self,dict):
        '''
            updates the attr attribute with values from the dictionary

            Parameters
            ----------
            dict : python dictionary
                Key value pairs to update

            Returns
            -------
            None
        '''
        self.attr.update(dict)
        return self


    def as_record(self):
        '''
            Returns the Locus as a record.
            NOTE: does not include Locus attributes. See Locus.as_dict()

            Parameters
            ----------
            None

            Returns
            -------
            A tuple containing Locus information
        '''
        return (self.chrom,self.start,self.end,self.name,self.window,self.global_order,self.feature_order,self.feature_type,self.id)

    @classmethod
    def from_record(cls,tpl):
        '''
            Creates a Locus object from a record.

            Parameters
            ----------
            record : tuple
                Tuple containing Locus information.

            Returns
            -------
            Locus object

        '''
        return cls(*tpl)

    def __setitem__(self,key,val):
        '''
            Set a Locus attribute.

            Parameters
            ----------
            key : str
                Attribute name.
            val : object
                Attribute value.

            Returns
            -------
            None
        '''
        self.attr[str(key)] = val

    def __getitem__(self,key):
        '''
            Retrieve a Locus attribute.

            Parameters
            ----------
            key : str
                Attribute name.

            Returns
            -------
            val : object
                Attribute value.
        '''
        return self.attr[key]

    def default_getitem(self,key,default=None):
        '''
            Return a default value if the attr[key] value is None

            Parameters
            ----------
            key : str
                Attribute name.
            default : str (default: None)
                A default value to return if the key is not present
                in the Locus attribute table.

            Returns
            -------
            val : object
                Attribute value of default value.
        '''
        if key not in self.attr.keys() or self[key] is None:
            return default
        else:
            return self[key]

    @property
    def start(self):
        '''
            Return the locus start value.

            Parameters
            ----------
            None

            Returns
            -------
            The Locus start position.
            NOTE: the minimum return value is 0
        '''
        return max(0,int(self._start))

    @property
    def end(self):
        '''
            Return the locus end value.

            Parameters
            ----------
            None

            Returns
            -------
            The locus end position.

        '''
        return int(self._end)

    @property
    def coor(self):
        '''
            Returns Locus start and stop positions

            Parameters
            ----------
            None

            Returns
            -------
            A len(2) tuple containing the start and stop positions

        '''
        return (self.start,self.end)
    @property
    def global_order(self):
        '''
            Returns the chromosome order of the locus regardless of feature type

            Parameters
            ----------
            None

            Returns
            -------
            The chromosome order of the locus in global scope.
        '''
        return int(self._global_order)

    @property
    def feature_order(self):
        '''
            Returns the chromosome order of the locus with respect to feature type

            Parameters
            ----------
            None

            Returns
            -------
            The chromosome order of the locus in narrow scope.
        '''
        return int(self._feature_order)

    @property
    def feature_type(self):
        '''
            Returns the feature type of the locus object

            Parameters
            ----------
            None

            Returns
            -------
            A string representing the feature type (gene, transposable element, snp, etc.)
        '''
        return str(self._feature_type)

    @property
    def upstream(self):
        '''
            Returns the starting position considered upstream of the Locus
            NOTE: Locus.window must be set for value to be different from
                  the locus start position.

            Parameters
            ----------
            None

            Returns
            -------
            The upstream position from the locus (based on locus window)

        '''
        return self.start - self.window

    @property
    def downstream(self):
        '''
            Returns the starting position considered downstream of the Locus
            NOTE: Locus.window must be set for value to be different from
                  the locus end position.

            Parameters
            ----------
            None

            Returns
            -------
            The downstream position from the locus (based on locus window)

        '''
        return self.end + self.window

    @property
    def name(self):
        '''
            A convenience method to access a Locus id

            Parameters
            ----------
            None

            Returns
            -------
            The locus id (name)
        '''
        return self.id

    @property
    def center(self):
        '''
            Returns the center position of the Locus.
            NOTE: this is a class property

            Parameters
            ----------
            None

            Returns
            -------
            A coordinate representing the center of the locus
        '''
        return (self.start + self.end)/2


    def center_distance(self, locus):
        '''
            Returns the distance between two Loci's center

            Parameters
            ----------
            locus: a object of type Locus

            Returns
            -------
            The distance between the loci centers
        '''
        if self.chrom == locus.chrom:
            return self.center - locus.center
        else:
            return np.inf


    def __add__(self,locus):
        '''
            Collapse two loci into a new 'meta' locus. The start is the min
            between the two loci and the window extends within as well as 1/2
            upstream and downstream of the original window sizes

            e.g. Locus(1,10,20) + Locus(1,30,40)

            Parameters
            ----------
            locus : another Locus object

            Returns
            -------
            a **new** Locus object containing the original loci as
            sub-loci

        '''
        # must be on the same chromosome to collapse
        if self-locus == float('Inf'):
            raise TypeError('Loci must be on same chromosome to collapse.')
        new_start = int(max(0,min(self.start,locus.start)))
        new_end = int(max(self.end,locus.end))
        new_window = self.window
        # This can be a list, since set gets called anyways
        new_sub_loci = self.sub_loci | locus.sub_loci

        return Locus(
            self.chrom, new_start, new_end,
            window=new_window, sub_loci=new_sub_loci
        )

    def __eq__(self,locus):
        '''
           Operator to compare equality of two Loci. They must share the
           same coordinates (chrom, start, end), everything else can differ.
           i.e.: Names or attrs can differ b/w loci

           e.g. Locus(1,10,100) == Locus(1,10,200) -> False
                Locus(1,10,100) == Locus(1,10,100) -> True

           Parameters
           ----------
           locus : another Locus object

           Returns
           -------
           bool : True if the same, else False
        '''
        if (self.chrom == locus.chrom and
            self.start == locus.start and
            self.end == locus.end):
            return True
        else:
            return False

    def within(self,locus):
        '''
            >>> x.within(y)

            Returns True if the coordinates of x are completely within
            the coordinates of y

            True:
            x:      |----===----|
            y:  |-------========--------|

            Parameters
            ----------
            locus: object of type Locus

            Returns
            -------
            bool
        '''
        if (locus.chrom == self.chrom
           and self.upstream >= locus.upstream
           and self.downstream <= locus.downstream):
            return True
        else:
            return False

    def encloses(self,locus):
        '''
            >>> x.encloses(y)

            Returns True if the start position of x is upstream
            and the end position of x is downstream of y.

            True:
            x:  |-------========--------|
            y:      |----===----|
        '''
        if (locus.chrom == self.chrom
           and locus.upstream >= self.upstream
           and locus.downstream <= self.downstream):
            return True
        else:
            return False

    def __contains__(self,locus):
        '''
            Test if a locus is within another locus, i.e. overlapping.
            NOTE: this function includes the locus **window** in the
                  calculation, the loci do not need to properly overlap,
                  however this case can be coverered by setting the
                  locus window to 0.

            Parameters
            ----------
            locus : another Locus object

            Returns
            -------
            bool : True or False

        '''
        if (locus.chrom == self.chrom and
               # The locus has as 'start' position within the Locus window
               # self:      |----===----|
               # locus:           |----===----|
               (( locus.upstream >= self.upstream and locus.upstream <= self.downstream)
               # The locus has an 'end' position within the Locus window
               or(locus.downstream <= self.downstream and locus.downstream >= self.upstream)
               # self:        |----===----|
               # locus:  |----===----|
               or(locus.upstream <= self.upstream and locus.downstream >= self.downstream)
               # self:        |----===----|
               # locus:  |----============----|
               or(locus.upstream >= self.upstream and locus.downstream <= self.downstream)
               # self:        |----============----|
               # locus:           |----===----|
            )):
            return True
        else:
            return False


    def __len__(self):
        '''
            Return the inclusive length of locus

            Returns
            -------
            The length of the locus (plus 1 bp)

        '''
        if self.start == self.end:
            return 1
        else:
            return self.end - self.start

    # XXX I checked the python doc and this method should never be called because
    # the rich comparison operators are defined (i.e. __eq__, __lt__, __gt__)
    def __cmp__(self,locus):
        '''
            Operator to compare the positions of two loci.

            Parameters
            ----------
            locus : a Locus object

            Returns
            -------
            a positive or negative number based on orientation,
            a.k.a. what is expected by python to do sorting
        '''
        if self.chrom == locus.chrom:
            return self.start - locus.start
        elif self.chrom > locus.chrom:
            return 1
        else:
            return -1

    def __lt__(self,locus):
        '''
            Operator to test if a locus is upstream from another locus

            Parameters
            ----------
            locus : a locus, which is tested to be upstream

            Returns
            -------
            bool
        '''
        if self.chrom == locus.chrom:
            if self.start == locus.start:
                return self.end < locus.end
            return self.start < locus.start
        else:
            return self.chrom < locus.chrom

    def __gt__(self,locus):
        '''
            Operator to test if a locus is downstream from another locus

            Parameters
            ----------
            locus : a locus, which is tested to be downstream

            Returns
            -------
            bool

        '''
        if self.chrom == locus.chrom:
            return self.start > locus.start
        else:
            return self.chrom > locus.chrom

    def __sub__(self,other):
        '''
            Return the inter-locus distance between two loci.

            Parameters
            ----------
            locus : a locus object

            Returns
            -------
            int : the distance between two loci
        '''
        if self.chrom != other.chrom:
            return float('Inf')
        if self == other:
            return 0
        else:
            # sort them
            a,b = sorted([self,other])
            return b.start - a.end

    def __str__(self):
        '''
            Return a string representation of the Locus

            Returns
            -------
            str
        '''
        return '''<{}>{}:{}-{}+{}({})'''.format(
            self._id, self.chrom,
            self.start, self.end,
            self.window, len(self.sub_loci)-1
        )

    def summary(self):
        '''
            Return summary information about a locus

            Parameters
            ----------
            None

            Returns
            -------
            A printable string with Locus info
        '''
        return '\n'.join([
            'ID: {}',
            'Chromosome: {}',
            'Start Position: {}',
            'End Position: {}',
            'Window Size: {}',
            'Feature Order: {}',
            'Feature Type: {}',
            'Additional attributes: {}',
            'Sub Loci: {}'
        ]).format(
            self._id, self.chrom,
            self.start, self.end,
            self.window, self.feature_order,
            self.feature_type, len(self.attr),
            len(self.sub_loci)
        )

    def __repr__(self):
        '''
           A convenience method for iPython
        '''
        return str(self)

    def __hash__(self):
        '''
            Convert the locus to a hash, uses md5.

            Parameters
            ----------
            None

            Returns
            -------
            str : md5 hash of locus

        '''
        digest = hashlib.md5(
            str.encode(str(self))
        ).hexdigest()
        return int(digest,base=16)
