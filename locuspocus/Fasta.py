from collections import defaultdict
import logging

import re
import numpy as np

from minus80 import Freezable
import reprlib
import pprint

from functools import lru_cache
from locuspocus import Chromosome


class Fasta(Freezable):
    '''
        A pythonic interface to a FASTA file. This interface
        allows convenient slicing into contigs (chromosomes).

       >>> from locuspocus import Fasta
       >>> x = Fasta.from_file('example.fa')

    '''
    log = logging.getLogger(__name__)                                               
    handler = logging.StreamHandler()                                               
    formatter = logging.Formatter(                                                  
                    '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'           
                )                                                                   

    handler.setFormatter(formatter)                                                 
    if not len(log.handlers): 
        log.addHandler(handler)                                                         
        log.setLevel(logging.INFO)          

    def __init__(self,name):
        '''
            Load a Fasta object from the Minus80.

            Parameters
            ----------
            name : str
                The name of the frozen object

            Returns
            -------
            A Fasta object
        '''
        super().__init__(name)
        # Load up from the database
        self._initialize_tables()

    def _initialize_tables(self):
        '''
            Initialize the tables for the FASTA class

            NOTE: internal method
        '''
        cur = self._db.cursor()

        cur.execute('''
            CREATE TABLE IF NOT EXISTS added_order (
                aorder INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT
            );
        ''')

        cur.execute('''
            CREATE TABLE IF NOT EXISTS nicknames (
                nickname TEXT,
                chrom TEXT,
                PRIMARY KEY(nickname,chrom),
                FOREIGN KEY(chrom) REFERENCES chroms(chrom)
            ) 
        ''')

        cur.execute('''
            CREATE TABLE IF NOT EXISTS attributes (
                chrom TEXT,
                attribute TEXT,
                PRIMARY KEY(chrom,attribute),
                FOREIGN KEY(chrom) REFERENCES chroms(chrom)
            )
        ''')

    def add_chrom(self,chrom,cur=None,force=False):
        '''
            Add a chromosome to the Fasta object.

            Parameters
            ----------
            name : str
                The name of the chromosome
        '''
        # Check for duplicates
        if chrom.name in self:
            if not force:
                raise ValueError(f'{chrom.name} already in FASTA')
        else:
            if cur is None:
                cur = self._db.cursor()
            cur.execute(
                '''
                INSERT OR REPLACE INTO added_order 
                    (name) 
                VALUES (?)
                ''',(chrom.name,)
            )
            for x in chrom._attrs:
                self._add_attribute(chrom.name,x)
        seqarray = np.array(chrom.seq)
        self.log.info(f'Adding {chrom.name}') 
        self._bcolz_array(chrom.name,seqarray)
        self.cache_clear()

    def cache_clear(self):
        self.__getitem__.cache_clear()

    @classmethod
    def from_file(cls,name,fasta_file,force=False):
        '''
            Create a Fasta object from a file.
        '''    
        self = cls(name)
        with open(fasta_file,'r') as IN, self._db as db: 
            cur = db.cursor()
            cur_chrom = None
            for line in IN:
                line = line.strip()
                if line.startswith('>'):
                    # Finish the last chromosome before adding a new one
                    if cur_chrom:
                        self.add_chrom(cur_chrom,cur=cur,force=force)
                    name,*attrs = line.lstrip('>').split()
                    cur_chrom = Chromosome(name,'',*attrs)
                else:
                    cur_chrom.seq = np.append(cur_chrom.seq,list(line))
            # Add the last chromosome
            self.add_chrom(cur_chrom,force=force,cur=cur)
        return self

    def __iter__(self):
        '''
            Iterate over chromosome objects
        '''
        chroms  = self._db.cursor().execute('SELECT name FROM added_order ORDER BY aorder')
        for (chrom,) in chroms:
            yield self[chrom]

    def __len__(self):
        '''
            Returns the number of chroms in the Fasta
        '''
        return self._db.cursor().execute('''
            SELECT COUNT(*) FROM added_order
        ''').fetchone()[0]

    def __contains__(self,obj):
        '''
           Returns boolean indicating if a named
           contig (chromosome) is in the fasta.
        '''
        if isinstance(obj,Chromosome):
            obj = obj.name 
        cur = self._db.cursor()
        # Check if in chrom names
        in_added = cur.execute(''' 
            SELECT COUNT(*) FROM added_order
            WHERE name = ?
        ''',(obj,)).fetchone()[0]
        if in_added == 1:
            return True
        # Check if in aliases
        in_alias = cur.execute('''
            SELECT COUNT(*) FROM nicknames
            WHERE nickname = ?
        ''',(obj,)).fetchone()[0]
        if in_alias == 1:
            return True
        # Otherise its not here
        return False
  

    @lru_cache(maxsize=128)
    def __getitem__(self,chrom_name):
       try:
            seq_array = self._bcolz_array(chrom_name)
       except Exception as e:
            chrom_name = self._get_nickname(chrom_name)
            seq_array = self._bcolz_array(chrom_name)
       finally:
            attrs = [x[0] for x in self._db.cursor().execute('''
                SELECT attribute FROM attributes 
                WHERE chrom = ?
            ''',(chrom_name,))]
            return Chromosome(chrom_name,seq_array,*attrs)


    def _add_attribute(self,chrom_name,attr,cur=None):
        '''
            Add an attribute the the Fasta object.
            Attributes describe chromosomes and 
            often follow the '>' token in the FASTA file.

            Parameters
            ----------
            chrom_name : str
                The name of the chromosome you are adding
                an attribute to
            attr : str
                the attribute you are adding
        '''
        if cur is None:
            cur = self._db.cursor()
        cur.execute(
            '''
            INSERT INTO attributes 
                (chrom,attribute) 
             VALUES (?,?)
             ''',
            (chrom_name,attr)
        )
        self.cache_clear()

    def _add_nickname(self,chrom,nickname,cur=None):
        '''
            Add a nickname for a chromosome

            Parameters
            ----------
            chrom : str
                The chromosome you want to nickname
            nickname : str
                The alternative name for the chromosome
        '''
        if cur is None:
            cur = self._db.cursor()
        cur.execute(
            '''
            INSERT OR REPLACE INTO nicknames 
                (nickname,chrom) 
            VALUES (?,?)
            ''',
            (nickname,chrom)
        )

    def _get_nickname(self,nickname):
        '''
            Get a chromosomem name by nickname
        '''
        return self._db.cursor().execute('''
            SELECT chrom FROM nicknames 
            WHERE nickname = ?
        ''',(nickname,)).fetchone()[0]

    def __repr__(self): #pragma: nocover
        return pprint.saferepr(
            reprlib.repr(list(self))
        )
