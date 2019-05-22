import pytest
import numpy as np

from itertools import chain
from locuspocus import Locus

from locuspocus.Exceptions import StrandError,ChromosomeError

@pytest.fixture
def simple_Locus():
    return Locus(1,100,200,attrs={'foo':'bar'}) 

def test_initialization(simple_Locus):
    # numeric chromosomes
    assert simple_Locus.chromosome == 1
    assert simple_Locus.start == 100
    assert simple_Locus.end == 200
    assert len(simple_Locus) == 101

def test_as_dict(simple_Locus):
    simple_Locus_dict = {
            'chromosome': 1, 
            'start': 100, 
            'end': 200, 
            'source': 'locuspocus', 
            'feature_type': 'locus', 
            'strand': '+', 
            'frame': None, 
            'name': None, 
            '_frozen': True, 
            '_LID': None
    }
    assert simple_Locus.as_dict() == simple_Locus_dict

def test_frozen(simple_Locus):
    try:
        # This should throw an exception
        simple_Locus.chromosome = 2
    except Exception as e:
        assert True

def test_setitem(simple_Locus):
    try:
        # This should throw an exception
        simple_Locus['foo'] = 'baz'
    except Exception as e:
        return True

def test_getitem(simple_Locus):
    assert simple_Locus['foo'] == 'bar'

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

def test_eq(simple_Locus):
    another_Locus = Locus(1, 110, 220)
    assert simple_Locus == simple_Locus
    assert simple_Locus != another_Locus

def test_loci_eq(testRefGen):
    '''
        Test the equals operator
    '''
    x = testRefGen['GRMZM2G018447']
    y = testRefGen['GRMZM2G018447']
    assert x == y

def test_loci_lt_by_chrom(testRefGen):
    x = Locus('1',1,1)
    y = Locus('2',1,1)
    assert x < y

def test_loci_gt_by_chrom(testRefGen):
    x = Locus('1',1,1)
    y = Locus('2',1,1)
    assert y > x

def test_loci_lt_by_pos(testRefGen):
    x = Locus('1',1,100)
    y = Locus('1',2,100)
    assert x < y

def test_loci_gt_by_pos(testRefGen):
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
    assert repr(simple_Locus) == "Locus(chromosome=1, start=100, end=200, source='locuspocus', feature_type='locus', strand='+', frame=None, name=None)"

def test_hash(simple_Locus):
    assert hash(simple_Locus) == 471250808972275261

def test_subloci_getitem_db_backed(testRefGen):
    x = testRefGen['GRMZM2G176595']
    assert x.subloci[0].name == 'GRMZM2G176595_T01'

def test_subloci_getitem():
    x = Locus('1',1,2)
    y = Locus('1',3,4,name='sublocus')
    x.add_sublocus(y)
    assert x.subloci[0].name == 'sublocus'

def test_subloci_iter():
    x = Locus('1',1,2)
    y = Locus('1',3,4,name='sublocus1')
    z = Locus('1',3,4,name='sublocus2')

    x.add_sublocus(y)
    x.add_sublocus(z)
    for sub in x.subloci:
        assert sub.chromosome == '1' 

def test_subloci_iter_db_backed(testRefGen):
    x = testRefGen['GRMZM2G008687']
    for sub in x.subloci:
        assert sub.feature_type == 'mRNA'

def test_subloci_len_db_backed(testRefGen):
    x = testRefGen['GRMZM2G008687']
    assert len(x.subloci) == 3

def test_subloci_len():
    x = Locus('1',1,2)
    y = Locus('1',3,4,name='sublocus1')
    z = Locus('1',3,4,name='sublocus2')

    x.add_sublocus(y)
    x.add_sublocus(z)
    assert len(x.subloci) == 2

def test_no_append_db_backed_subloci(testRefGen):
    x = testRefGen['GRMZM2G008687']
    with pytest.raises(ValueError):
        x.subloci.append(0)

def test_attrs_dbgetitem_on_non_db_backed(testRefGen):
    x = Locus('1',3,4,attrs={'foo':'locus1'})
    with pytest.raises(ValueError):
        x.attrs._db_getitem('foo')

def test_attrs_dbgetitem(testRefGen):
    x = testRefGen['GRMZM2G008687']
    assert x.attrs['Name'] == 'GRMZM2G008687' 

def test_attrs_dbgetitem_not_there(testRefGen):
    x = testRefGen['GRMZM2G008687']
    with pytest.raises(KeyError):
        assert x.attrs['NO'] 

def test_attrs_cannot_set_item(testRefGen):
    x = testRefGen['GRMZM2G008687']
    with pytest.raises(ValueError):
        x.attrs['foo'] = 'bar'

def test_attrs_keys_method_db_backed(testRefGen):
    x = testRefGen['GRMZM2G008687']
    assert sorted(x.attrs.keys()) == ['ID', 'Name', 'biotype']

def test_attrs_keys_method(testRefGen):
    x = Locus('1',3,4,attrs={'foo':'locus1','bar':'baz'})
    assert sorted(x.attrs.keys()) == ['bar','foo']

def test_attrs_keys_method_empty(testRefGen):
    x = Locus('1',3,4,attrs={})
    assert len(list(x.attrs.keys())) == 0

def test_attrs_vals_method_db_backed(testRefGen):
    x = testRefGen['GRMZM2G008687']
    assert len(sorted(x.attrs.values())) == 3

def test_attrs_vals_method(testRefGen):
    x = Locus('1',3,4,attrs={'foo':'locus1','bar':'baz'})
    assert len(sorted(x.attrs.values())) == 2

def test_attrs_vals_method_empty(testRefGen):
    x = Locus('1',3,4,attrs={})
    assert len(list(x.attrs.values())) == 0

def test_attrs_items_db_backed(testRefGen):
    x = testRefGen['GRMZM2G008687']
    assert len(sorted(x.attrs.items())) == 3

def test_attrs_items(testRefGen):
    x = Locus('1',3,4,attrs={'foo':'locus1','bar':'baz'})
    assert len(sorted(x.attrs.items())) == 2

def test_attrs_contains(testRefGen):
    x = Locus('1',3,4,attrs={'foo':'locus1','bar':'baz'})
    assert 'foo' in x.attrs

def test_attrs_contains_db_backed(testRefGen):
    x = testRefGen['GRMZM2G008687']
    assert 'biotype' in x.attrs

def test_attrs_contains_db_backed_missing(testRefGen):
    x = testRefGen['GRMZM2G008687']
    assert 'abcdefg' not in x.attrs

def test_attrs_repr():
    x = Locus('1',3,4,attrs={'foo':'locus1','bar':'baz'})
    assert repr(x.attrs)

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

def test_as_record():
    x = Locus('1',3,4,strand='+')
    # This doesn't compare the dictionaries of each ...
    assert x.as_record()[0] == ('1', 3, 4, 'locuspocus', 'locus', '+', None, None, 1030788771219373542)

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
    assert x in z.subloci
    assert y in z.subloci

def test_combine_chromosome_mismatch():
    x = Locus('1',1,2)
    y = Locus('2',3,4)
    with pytest.raises(ChromosomeError):
        z = x.combine(y)
