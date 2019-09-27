import pytest
from locuspocus import Locus,Loci,Term
from locuspocus.exceptions import *
import minus80 as m80
from dataclasses import FrozenInstanceError

@pytest.fixture
def t():
    a = Locus('1',1,100)
    b = Locus('2',200,300)
    c = Locus('3',300,400)
    t = Term(
        'test',
        desc='test term',
        loci = [a,b,c],
        attrs = {
            'foo' : 'bar',
            'baz' : 'bat'
        }
    )
    return t

def test_term_init(t):
    assert t

def test_term_name(t):
    assert t.name == 'test'

def test_term_desc(t):
    assert t.desc == 'test term'

def test_term_loci(t):
    for l in t.loci:
        assert isinstance(l,Locus)
        
def test_term_attrs(t):
    assert t.attrs['foo'] == 'bar'
    assert t.attrs['baz'] == 'bat'

def test_no_parent(t):
    assert t.parent is None

def test_frozen(t):
    assert t._frozen == True

def test_len(t):
    assert len(t) == 3

def test_getitem(t):
    assert t['foo'] == 'bar'
    assert t['baz'] == 'bat'

def test_setitem_frozen(t):
    with pytest.raises(FrozenInstanceError):
        # throws an error since Terms are frozen
        t['fizz'] = 'wizz'

def test_add_locus():
    t = Term('test')
    assert len(t) == 0
    t.add_locus(Locus('1',1,100))
    assert len(t) == 1

def test_add_duplicate_locus():
    t = Term('test')
    assert len(t) == 0
    t.add_locus(Locus('1',1,100))
    assert len(t) == 1
    # Add duplicate which wont count
    t.add_locus(Locus('1',1,100))
    assert len(t) == 1

def test_nearby_loci(t):
    x = Locus('1',150,200)
    nearby = list(t.nearby_loci(x,max_distance=50))
    assert len(nearby) == 1
    nearby = list(t.nearby_loci(x,max_distance=10))
    assert len(nearby) == 0
    
