from typing import (
    Generator
)

class SubLoci():
    # A restricted list interface to subloci
    def __init__(self,loci=None):
        self._loci = loci

    @property
    def empty(self) -> bool:
        if self._loci is None:
            return True
        else:
            return False

    def __eq__(self, other) -> bool:
        if len(self) != len(other):
            return False
        else:
            return sorted(self) == sorted(other)

    def __iter__(self) -> Generator['Locus',None,None]:
        if self.empty:
            return (x for x in [])
        return (x for x in self._loci)

    def add(self, locus: 'Locus') -> None:
        if self.empty:
            self._loci = []
        self._loci.append(locus)

    def __getitem__(self,index: int) -> "Locus":
        if self.empty:
            raise IndexError
        return self._loci[index]

    def __len__(self) -> int:
        if self.empty:
            return 0
        return len(self._loci)

    def __repr__(self) -> str:
        if len(self) == 0 or self.empty:
            return "[]"
        return '[\n'+'\t\n'.join([repr(x) for x in self])+'\n]'

    def find(
        self,
        name
    ) -> "Locus":
        '''
            Find a sublocus by its name.
        '''
        if self.empty:
            return None
        for l in self.traverse():
            if l.name == name:
                return l
        return None

    def traverse(
        self,
        mode='depth'
    ) -> Generator['Locus',None,None]:
        '''
            Perform a depth first traversal of subloci
        '''
        if self.empty:
            return None
        for l in self:
            yield l
            yield from l.subloci.traverse()


