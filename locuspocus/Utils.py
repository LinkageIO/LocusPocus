import re

def guess_type(object):
    '''
        Guess the type of object from the class attribute
    '''
    # retrieve a list of classes
    classes = re.match(
        "<class '(.+)'>",
        str(object.__class__)
    ).groups()[0].split('.')
    # Return the most specific one
    return classes[-1]
