import logging
import re
import reprlib
import pprint

import numpy as np

from minus80 import Freezable
from minus80.RawFile import RawFile
from collections import defaultdict
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

    def __init__(self,name,parent=None):
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
        super().__init__(name,parent=parent)
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

    def add_chrom(self,chrom,replace=False,cur=None):
        '''
            Add a chromosome to the Fasta object.

            Parameters
            ----------
            chrom : Chromosome object
                The chromosome object to add.
                See LocusPocus.Chromosome
            replace : bool (default: False)
                By default a chromosome can only be added 
                once. If this is set, the chromosome object
                will be replaced.

        '''
        self.log.info(f'Adding {chrom.name}') 
        # Check for duplicates
        if chrom.name in self:
            if not replace:
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
        self._bcolz_array(chrom.name,seqarray)
        self.cache_clear()

    def del_chrom(self,chrom):
        '''
            Delete a chromosome from the database
        '''
        if isinstance(chrom,Chromosome):
            name = chrom.name
        elif isinstance(chrom,str):
            name = chrom
        else:
            raise ValueError(f'input must be a Chromosome object or a string')
        if name not in self:
            raise ValueError(f"'{name}' not in the {self._m80_dtype}('{self._m80_name}')")
        self._db.cursor().execute('''
            DELETE FROM added_order WHERE name = ?;
            DELETE FROM nicknames WHERE chrom = ?;
            DELETE FROM attributes WHERE chrom = ?;
        ''',(name,name,name))
        self._bcolz_remove(name)
        

    def chrom_names(self):
        '''
            Returns an iterable of chromosome names

            Parameters
            ----------
            None

            Returns
            -------
            An iterable of chromosome names in added order
        '''
        return (x for (x,) in self._db.cursor().execute('''
            SELECT name FROM added_order ORDER BY aorder
        '''))

    def cache_clear(self):
        self.__getitem__.cache_clear()

    @classmethod
    def from_file(cls,name,fasta_file,replace=False,parent=None):
        '''
            Create a Fasta object from a file.
        '''    
        self = cls(name,parent=parent)
        with RawFile(fasta_file) as IN, self._db as db: 
            cur = db.cursor()
            cur_chrom = None
            seqs = []
            name, attrs = None,None
            for line in IN:
                line = line.strip()
                if line.startswith('>'):
                    # Finish the last chromosome before adding a new one
                    if len(seqs) > 0:
                        cur_chrom = Chromosome(name,seqs,*attrs)
                        self.add_chrom(cur_chrom,cur=cur,replace=replace)
                        seqs = []
                    name,*attrs = line.lstrip('>').split()
                else:
                    seqs += line
                    #cur_chrom.seq = np.append(cur_chrom.seq,list(line))
            # Add the last chromosome
            cur_chrom = Chromosome(name,seqs,*attrs)
            self.add_chrom(cur_chrom,cur=cur,replace=replace)
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
        if chrom_name not in self:
            raise ValueError(f'{chrom_name} not in {self._m80_name}')
        try:
            seq_array = self._bcolz_array(chrom_name)
        except Exception as e:
            chrom_name = self._get_nickname(chrom_name)
            seq_array = self._bcolz_array(chrom_name)
        finally:
            attrs = [x[0] for x in self._db.cursor().execute('''
                SELECT attribute FROM attributes 
                WHERE chrom = ?
                ORDER BY rowid -- This preserves the ordering of attrs
            ''',(chrom_name,))]
            return Chromosome(chrom_name,seq_array,*attrs)

    def to_fasta(self,filename,line_length=70):
        '''
            Print the chromosomes to a file in FASTA format

            Paramaters
            ----------
            filename : str
                The output filename
            line_length : int (default: 70)
                The number of nucleotides per line

            Returns
            -------
            None
        '''
        with open(filename,'w') as OUT:                                                                                                        
            for chrom_name in self.chrom_names():                                        
                print(f'Printing out {chrom_name}')
                chrom = self[chrom_name]                                                 
                #easy_id = ids[chrom_name]                                               
                start_length = len(chrom)
                #if easy_id == 'chrUn':
                #    easy_id = easy_id + '_' + chrom_name
                print(f'>{chrom_name} {"|".join(chrom._attrs)}',file=OUT)     
                printed_length = 0
                for i in range(0,len(chrom),70):                                        
                    sequence = chrom.seq[i:i+70]
                    print(''.join(sequence),file=OUT) 
                    printed_length += len(sequence)
                if printed_length != start_length: #pragma: no cover
                    raise ValueError('Chromosome was truncated during printing')
        return None

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
