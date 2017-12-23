import reprlib
import numpy as np

class Chromosome(object) :                                                          
    '''                                                                             
    A Chromosome is a lightweight object which maps indices to                     
    string positions. It's a named sequence. 

    NOTE: chromosome indices are different that python indices.
          Mainly, they are 1 indexed AND they are inclusive 
          for start and stop positions.

    >>> from locuspocus import Chromosome
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
    def __init__(self,name,seq,*args):
        self.name = name
        if isinstance(seq,str):
            self.seq = np.array(list(seq))
        else:
            self.seq = np.array(seq)
        self._attrs = set(args)


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
            reprlib.repr(''.join(self.seq[1:100]))
        )

    def __eq__(self,obj):
        if self.name == obj.name and all(self.seq == obj.seq):
            return True
        else:
            return False


