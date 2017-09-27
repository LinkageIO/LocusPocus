from collections import defaultdict
import logging

import re

from minus80 import Freezable

class Chromosome(object) :                                                          
    '''                                                                             
    A Chromosome is a lightweight object which maps indices to                     
    string positions.

    NOTE: chromosome indices are different that python indices.
          Mainly, they are 1 indexed AND they are inclusive 
          for start and stop positions.

    In [0]: from locuspocus.Fasta import Chromosome

    In [1]: x = Chromosome('AAACCCTTTGGG')

    In [3]: [1]
    Out[3]: [1]

    In [4]: x[1]
    Out[4]: 'A'

    In [5]: x[1:5]
    Out[5]: 'AAACC'

    In [6]: x[5:10]
    Out[6]: 'CCTTTG'

    In [7]: len(x)
    Out[7]: 12

    '''                                                                            
    def __init__(self,seq):                                                        
        self.seq = str(seq)                                                        
    def __getitem__(self,pos):                                                     
        if isinstance(pos,slice):                                                  
            if pos.start < 1:
                raise ValueError('Genetic coordinates cannot start less than 1')
            return self.seq[max(0,pos.start-1):pos.stop]                                  
        # chromosomes start at 1, python strings start at 0                         
        else:
            return self.seq[int(pos)-1]                                                
    def __len__(self):                                                             
        return len(self.seq)        

class Fasta(object):
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

    def __init__(self):
        '''
            Initialize an empty Fasta object.
        '''
        self._file = None       
        self.added_order = []
        self.chroms = {}
        self.nicknames = {}
        self.attributes = defaultdict(list)
        self._chrom_set = set(self.chroms.keys())


    def to_minus80(self,name):
        '''
            Store the Fasta Object in the Minus80

            Parameters
            ----------
            name : str
                The name of the store in the  minus80.
        '''
        from minus80 import Freezable
        store = Freezable(name,type=Freezable.guess_type(self))
        self._initialize_tables(store)
        cur = store._db.cursor()
        # Store the file
        store._dict('_file',self._file)
        # Store the chroms in added order
        for i,chrom in enumerate(self.added_order):
            seq = self.chroms[chrom].seq
            cur.execute('INSERT OR REPLACE INTO chroms (added_order,chrom,seq) VALUES (?,?,?)',(i,chrom,seq))
        for nickname,chrom in self.nicknames.items():
            cur.execute('INSERT OR REPLACE INTO nicknames (nickname,chrom) VALUES (?,?)',(nickname,chrom))
        for chrom,attrs in self.attributes.items():
            for attr in attrs:
                cur.execute('INSERT OR REPLACE INTO attributes (chrom,attribute) VALUES (?,?)',(chrom,attr))

    @classmethod
    def from_minus80(cls,name):
        self = cls()
        from minus80 import Freezable
        store = Freezable(name,type=Freezable.guess_type(self))
        self._file = store._dict('_file')
        cur = store._db.cursor()
        for chrom,seq in cur.execute('SELECT chrom,seq FROM chroms ORDER BY added_order ASC'):
            self.add_chrom(chrom,Chromosome(seq))
        for nickname,chrom in cur.execute('SELECT nickname,chrom FROM nicknames'):
            self.nicknames[nickname] = chrom
        for chrom,attr in cur.execute('SELECT chrom,attribute FROM attributes'):
            self.add_attribute(chrom,attr)
        return self
        


    def __contains__(self,item):
        '''
            Returns boolean indicating if a named
            contig (chromosome) is in the fasta.
        '''
        if item in self._chrom_set:
            return True
        elif item in self.chroms:
            return True
        else:
            return False

    def __getitem__(self,chrom_name):
        try:
            if chrom_name in self.chroms:
                return self.chroms[chrom_name]
            if chrom_name in self.nicknames:
                return self.chroms[self.nicknames[chrom_name]]
        except Exception as e:
            raise KeyError(
                '{} not in Fasta: {}'.format(chrom_name,self._file)
            )

    def add_chrom(self,chrom_name,chromosome):
        '''
            Add a chromosome to the Fasta object.
        '''
        self.added_order.append(chrom_name)
        self.chroms[chrom_name] = chromosome

    def add_attribute(self,chrom_name,attr):
        '''
            Add an attribute the the Fasta object.
            Attributes describe chromosomes and 
            often follow the '>' token in the FASTA file.
        '''
        self.attributes[chrom_name].append(attr)


    @classmethod
    def from_file(cls,fasta_file,nickname=None):
        '''
            Create a Fasta object from a file.
        '''    
        self = cls()
        self._file = fasta_file
        with open(self._file,'r') as IN: 
            cur_chrom = None
            cur_seqs = []
            for line in IN:
                line = line.strip()
                if line.startswith('>'):
                    # Finish the last chromosome before adding a new one
                    if cur_chrom:
                        self.add_chrom(cur_chrom,Chromosome("".join(cur_seqs)))
                    name,*attrs = line.lstrip('>').split()
                    cur_chrom = name
                    cur_seqs = []
                    for attr in attrs:
                        self.add_attribute(cur_chrom,attr)
                    if nickname != None:
                        pattern,replace = nickname
                        alt = re.sub(pattern,replace,line)
                        if alt != line:
                            log('Found a nickname: mapping {} -> {}',alt,name)
                            self.nicknames[alt] = name
                            self.nicknames[name] = alt
                else:
                    cur_seqs.append(line)
            # Add the last chromosome
            self.add_chrom(cur_chrom,Chromosome("".join(cur_seqs)))
        return self


    @staticmethod
    def _initialize_tables(store):
        cur = store._db.cursor()

        cur.execute('''
            CREATE TABLE IF NOT EXISTS chroms (
                added_order INT,
                chrom TEXT,
                seq BLOB,
                PRIMARY KEY(chrom)
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
