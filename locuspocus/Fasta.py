from collections import defaultdict
import logging

import re
import numpy as np

from minus80 import Freezable
import reprlib
import pprint

from functools import lru_cache

class Chromosome(object) :                                                          
    '''                                                                             
    A Chromosome is a lightweight object which maps indices to                     
    string positions.

    NOTE: chromosome indices are different that python indices.
          Mainly, they are 1 indexed AND they are inclusive 
          for start and stop positions.

    >>> from locuspocus.Fasta import Chromosome
    >>> x = Chromosome('AAACCCTTTGGG')
    >>> x[1]
    'A'
    >>> x[1:5]
    'AAACC'
    >>> x[5:10]
    'CCTTTG'
    >>> len(x)
    12

    '''                                                                            
    def __init__(self,name,seq):                                                        
        self.name = name
        self.seq = seq                                                      
    def __getitem__(self,pos):                                                     
        if isinstance(pos,slice):                                                  
            if pos.start < 1:
                raise ValueError('Genetic coordinates cannot start less than 1')
            return self.seq[max(0,pos.start-1):pos.stop]                                  
        # chromosomes start at 1, python strings start at 0                         
        else:
            if pos < 1:
                raise ValueError('Genetic coordinates cannot start less than 1')
            return self.seq[int(pos)-1]                                                
    def __len__(self):                                                             
        return len(self.seq)         

    def __repr__(self):
        return 'Chromosome({})'.format(
            reprlib.repr(self.seq)
        )

    def __eq__(self,obj):
        if self.name == obj.name and self.seq == obj.seq:
            return True
        else:
            return False


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

    def add_chrom(self,name,chrom,cur=None):
        '''
            Add a chromosome to the Fasta object.

            Parameters
            ----------
            name : str
                The name of the chromosome
        '''
        if cur is None:
            cur = self._db.cursor()
        cur.execute(
            '''
            INSERT OR REPLACE INTO added_order 
                (name) 
            VALUES (?)
            ''',(name,)
        )
        seqarray = np.array(list(chrom.seq))
        self.log.warn(f'Adding {name}') 
        self._bcolz_array(name,seqarray)
        

    @classmethod
    def from_file(cls,name,fasta_file,nickname=None):
        '''
            Create a Fasta object from a file.
        '''    
        self = cls(name)
        with open(fasta_file,'r') as IN, self._db as db: 
            cur = db.cursor()
            cur_chrom = None
            cur_seqs = []
            for line in IN:
                line = line.strip()
                if line.startswith('>'):
                    # Finish the last chromosome before adding a new one
                    if cur_chrom:
                        self.add_chrom(cur_chrom,Chromosome("".join(cur_seqs)),cur=cur)
                    name,*attrs = line.lstrip('>').split()
                    cur_chrom = name
                    cur_seqs = []
                    for attr in attrs:
                        self.add_attribute(cur_chrom,attr,cur=cur)
                    if nickname != None:
                        pattern,replace = nickname
                        alt = re.sub(pattern,replace,line)
                        if alt != line:
                            self.log('Found a nickname: mapping {} -> {}',alt,name)
                            self._add_nickname(name,alt,cur=cur)
                            self._add_nickname(alt,name,cur=cur)
                else:
                    cur_seqs.append(line)
            # Add the last chromosome
            self.add_chrom(cur_chrom,Chromosome("".join(cur_seqs)))
        return self

    def __iter__(self):
        '''
            Iterate over chromosome objects
        '''
        chroms  = self._db.cursor().execute('SELECT name FROM added_order ORDER BY aorder')
        for (chrom,) in chroms:
            yield self[chrom]

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
           return Chromosome(chrom_name,seq_array)
       except Exception as e:
           pass
       try:
           chrom_name = self._get_nickname(chrom_name)
           return self[chrom_name]
       except Exception as e:
           raise KeyError(
               '{} not in Fasta: {}'.format(chrom_name,self._m80_name)
           )

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
            INSERT OR REPLACE INTO attributes 
                (chrom,attribute) 
             VALUES (?,?)
             ''',
            (chrom_name,attr)
        )

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

    def __repr__(self):
        return pprint.saferepr(
            reprlib.repr(self.chroms)
        )
