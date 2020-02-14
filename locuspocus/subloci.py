

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

    def __iter__(self):
        if self.empty:
            return (x for x in [])
        return (x for x in self._loci)

    def add(self,locus) -> None:
        if self.empty:
            self._loci = []
        self._loci.append(locus)

    def __getitem__(self,index) -> "Locus":
        if self.empty:
            raise IndexError
        return self._loci[index]

    def __len__(self) -> int:
        if self.empty:
            return 0
        return len(self._loci)

    def __repr__(self):
        if len(self) == 0 or self.empty:
            return "None"
        return f"[{len(self)} subloci]"

    def find(
        self,
        name
    ) -> "Locus":
        '''
            Find a sublocus by its name.
        '''
        if self.empty:
            return None
        for l in self:
            if l.name == name:
                return l
            else:
                # Recurse
                if l.subloci.find(name) is None:
                    continue
        return None


