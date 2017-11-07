# Exception abstract class

class ZeroWindowError(Exception): # pragma: no cover
    def __init__(self,expr,message,*args):
        self.expr = expr
        self.message = (
            'Operation requiring window, but window is 0:' + \
            message.format(args)
        )
