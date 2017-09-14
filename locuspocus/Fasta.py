from collections import defaultdict
from ponytools.Tools import log
from .Chromosome import Chromosome

from .Exceptions import MissingChromosomeError

import re

class Fasta(object):
    def __init__(self):
        self._file = None
        self.added_order = []
        self.chroms = {}
        self.nicknames = {}
        self.attributes = defaultdict(list)
        self._chrom_set = set(self.chroms.keys())

    def __contains__(self,item):
        if item in self._chrom_set:
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
            raise MissingChromosomeError(
                '{} not in Fasta: {}'.format(chrom_name,self._file)
            )

    def add_chrom(self,chrom_name,chromosome):
        self.added_order.append(chrom_name)
        self.chroms[chrom_name] = chromosome

    def add_attribute(self,chrom_name,attr):
        self.attributes[chrom_name].append(attr)


    @classmethod
    def from_file(cls,fasta_file,nickname=None):
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



