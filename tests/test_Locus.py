import pytest
import numpy as np

from itertools import chain
from locuspocus import Locus

from locuspocus.Exceptions import StrandError,ChromosomeError,LocusError

@pytest.fixture
def simple_Locus():
    return Locus('1',100,200,attrs={'foo':'bar'}) 

def test_initialization(simple_Locus):
    # numeric chromosomes
    assert simple_Locus.chromosome == '1'
    assert simple_Locus.start == 100
    assert simple_Locus.end == 200
    assert len(simple_Locus) == 101

def test_bad_core_attribute(simple_Locus):
    with pytest.raises(AttributeError):
        simple_Locus.foo

def test_none_core_attribute_returns_none(simple_Locus):
    assert simple_Locus.frame is None

def test_getitem(simple_Locus):
    assert simple_Locus['foo'] == 'bar'

def test_getitem_int():
    x = Locus('1',100,200,attrs={'foo':12}) 
    assert isinstance(x['foo'],int)
def test_getitem_float():
    x = Locus('1',100,200,attrs={'foo':12.0}) 
    assert isinstance(x['foo'],float)
def test_getitem_str():
    x = Locus('1',100,200,attrs={'foo':'12'}) 
    assert isinstance(x['foo'],str)
def test_getitem_bool():
    x = Locus('1',100,200,attrs={'foo':True}) 
    assert isinstance(x['foo'],bool)

def test_bad_key_set_item():
    with pytest.raises(ValueError):
        x = Locus('1',100,200,attrs={1:True}) 

def test_bad_val_set_item():
    with pytest.raises(ValueError):
        x = Locus('1',100,200,attrs={'foo':{}}) 

def test_default_getitem(simple_Locus):
    assert simple_Locus.default_getitem('name', 'default') == 'default'

def test_start(simple_Locus):
    assert simple_Locus.start == 100

def test_plus_stranded_start():
    l = Locus('1',1,100,strand='+')
    assert l.stranded_start == 1 

def test_minus_stranded_start():
    l = Locus('1',1,100,strand='-')
    assert l.stranded_start == 100 

def test_end(simple_Locus):
    assert simple_Locus.end == 200

def test_plus_stranded_end():
    l = Locus('1',1,100,strand='+')
    assert l.stranded_end == 100 

def test_minus_stranded_end():
    l = Locus('1',1,100,strand='-')
    assert l.stranded_end == 1 

def test_coor(simple_Locus):
    assert simple_Locus.coor == (100, 200)

def test_upstream(simple_Locus):
    assert simple_Locus.upstream(50) == 50

def test_upstream_minus_strand(simple_Locus):
    l = Locus('1',1,100,strand='-')
    assert l.upstream(50) == 150

def test_downstream(simple_Locus):
    assert simple_Locus.downstream(50) == 250

def test_downstream_minus_strand(simple_Locus):
    l = Locus('1',100,200,strand='-')
    assert l.downstream(50) == 50

def test_center():
    l = Locus('1',100,200,strand='-')
    assert l.center == 150.5 

def test_name(simple_Locus):
    assert simple_Locus.name is None

def test_eq():
    x = Locus(1, 110, 220)
    y = Locus(1, 110, 220)
    assert x == y

def test_ne_diff_attrs():
    x = Locus(1, 110, 220,attrs={'foo':'bar'})
    y = Locus(1, 110, 220,attrs={'baz':'bat'})
    assert x != y

def test_loci_lt_by_chrom():
    x = Locus('1',1,1)
    y = Locus('2',1,1)
    assert x < y

def test_loci_gt_by_chrom():
    x = Locus('1',1,1)
    y = Locus('2',1,1)
    assert y > x

def test_loci_lt_by_pos():
    x = Locus('1',1,100)
    y = Locus('1',2,100)
    assert x < y

def test_loci_gt_by_pos():
    x = Locus('1',1,100)
    y = Locus('1',2,200)
    assert y > x

def test_len(simple_Locus):
    assert len(simple_Locus) == 101
    assert len(Locus(1, 100, 100)) == 1

def test_lt(simple_Locus):
    same_chrom_Locus = Locus(1, 110, 220)
    diff_chrom_Locus = Locus(2, 90, 150)
    assert simple_Locus < same_chrom_Locus
    assert simple_Locus < diff_chrom_Locus

def test_gt(simple_Locus):
    same_chrom_Locus = Locus(1, 90, 150)
    diff_chrom_Locus = Locus(2, 90, 150)
    assert simple_Locus > same_chrom_Locus
    assert diff_chrom_Locus > simple_Locus

def test_repr(simple_Locus):
    assert repr(simple_Locus) 

def test_repr_with_name():
    x = Locus('1',1,100,name='test')
    assert repr(x)

def test_children_len():
    x = Locus('1',1,2)
    y = Locus('1',3,4,name='sublocus1')
    z = Locus('1',3,4,name='sublocus2')

    x.children = (y,z)
    assert len(x.children) == 2

def test_children_detach_in_children_setter():
    x = Locus('1',1,2)
    y = Locus('1',3,4,name='sublocus1')
    z = Locus('1',3,4,name='sublocus2')

    x.children = [y]
    # Now trigger a detach for current children 
    # when assigning new children
    x.children = [z]
    assert len(x.children) == 1 


def test_attrs_keys_method_empty():
    x = Locus('1',3,4,attrs={})
    assert len(list(x.attrs.keys())) == 0

def test_attrs_vals_method_empty():
    x = Locus('1',3,4,attrs={})
    assert len(list(x.attrs.values())) == 0

def test_setitem():
    x = Locus('1',3,4,attrs={'foo':'locus1','bar':'baz'})
    assert x['foo'] == 'locus1'
    x['foo'] = 'bar'
    assert x['foo'] == 'bar'

def test_attrs_contains():
    x = Locus('1',3,4,attrs={'foo':'locus1','bar':'baz'})
    assert 'foo' in x.attrs

def test_le_equals():
    x = Locus('1',3,4,attrs={'foo':'locus1','bar':'baz'})
    y = Locus('1',3,4)
    assert x <= y

def test_le_less():
    x = Locus('1',3,4,attrs={'foo':'locus1','bar':'baz'})
    y = Locus('1',30,40)
    assert x <= y

def test_ge_equals():
    x = Locus('1',3,4,attrs={'foo':'locus1','bar':'baz'})
    y = Locus('1',3,4)
    assert x >= y

def test_ge_greater():
    x = Locus('1',30,40,attrs={'foo':'locus1','bar':'baz'})
    y = Locus('1',3,4)
    assert x >= y

def test_stranded_start_invalid():
    # Strand cannot be '='
    x = Locus('1',3,4,strand='=')
    with pytest.raises(StrandError):
        x.stranded_start

def test_stranded_stop_invalid():
    # Strand cannot be '='
    x = Locus('1',3,4,strand='=')
    with pytest.raises(StrandError):
        x.stranded_end

def test_center_distance():
    x = Locus('1',1,100,strand='+')
    # This needs to be 201 since x starts at 1
    y = Locus('1',201,300,strand='=')
    assert x.center_distance(y) == 200

def test_center_distance_different_chroms():
    x = Locus('1',1,100,strand='+')
    # This needs to be 201 since x starts at 1
    y = Locus('2',201,300,strand='+')
    assert x.center_distance(y) == np.inf

def test_str():
    x = Locus('1',1,100,strand='+')
    assert str(x) == repr(x)

def test_combine():
    x = Locus('1',1,2)
    y = Locus('1',3,4)
    z = x.combine(y)

    assert z.start == 1
    assert z.end == 4
    assert len(z.children) == 2

def test_combine_chromosome_mismatch():
    x = Locus('1',1,2)
    y = Locus('2',3,4)
    with pytest.raises(ChromosomeError):
        z = x.combine(y)

def test_distance():
    x = Locus('1',1,100)
    y = Locus('1',150,250)
    assert x.distance(y) == 49

def test_distance_diff_chroms():
    x = Locus('1',1,100)
    y = Locus('2',150,250)
    assert x.distance(y) == np.inf


def test_get_parent():
    x = Locus('1',1,2)
    y = Locus('1',3,4)
    
    x.parent = y
    assert x.parent == y

def test_get_none_parent():
    x = Locus('1',1,2)
    assert x.parent is None


def test_set_parent_bad_type():
    x = Locus('1',1,2)
    with pytest.raises(LocusError):
        x.parent = 12


