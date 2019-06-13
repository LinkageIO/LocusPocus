# Exception abstract class

class Error(Exception):
    pass

class ZeroWindowError(Exception): # pragma: no cover
    def __init__(self,expr,message,*args):
        self.expr = expr
        self.message = (
            'Operation requiring window, but window is 0:' + \
            message.format(args)
        )

class LocusError(Exception): # pragma: no cover
    def __init__(self,expr=None,message=None):
        self.expr = expr
        self.message = message

class MissingLocusError(Exception): # pragma: no cover
    def __init__(self,expr=None,message=None):
        self.expr = expr
        self.message = (
            f'Locus not present in RefLoci:{message}'
        )

class StrandError(Exception): #pragma: no cover
    def __init__(self,expr=None,message=''):
        self.expr = expr
        self.message = (
            f'Locus has a bad strand:{message}'
        )

class ChromosomeError(Exception): #pragma: no cover
    def __init__(self,message=''):
        self.message = (
            f'{message}'
        )

