class LocusAttrs:
    # a restricted dict interface to attributes
    def __init__(self, attrs=None):
        self._attrs = attrs

    def __len__(self):
        if self.empty:
            return 0
        else:
            return len(self._attrs)

    def __eq__(self, other):
        if self.empty and other.empty:
            return True
        elif len(self) != len(other):
            # Short circuit on length
            return False
        else:
            return sorted(self.items()) == sorted(other.items())

    @property
    def empty(self):
        if self._attrs is None:
            return True
        else:
            return False

    def keys(self):
        if self.empty:
            return []
        else:
            return self._attrs.keys()

    def values(self):
        if self.empty:
            return []
        else:
            return self._attrs.values()

    def items(self):
        if self.empty:
            return {}
        else:
            return self._attrs.items()

    def __contains__(self, key):
        if self.empty:
            return False
        return key in self._attrs

    def __getitem__(self, key):
        if self.empty:
            raise KeyError()
        return self._attrs[key]

    def __setitem__(self, key, val):
        if self.empty:
            self._attrs = {}
        self._attrs[key] = val

    def __repr__(self):
        if self.empty:
            return repr({})
        return repr(self._attrs)
