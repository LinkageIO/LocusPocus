import pytest
from locuspocus import Chromosome


@pytest.fixture
def chr1():
    return Chromosome('chr1','A'*500000)                                                  

@pytest.fixture
def chr2():
    return Chromosome('chr1','C'*500000)                                                  


def test_init(chr1):
    assert chr1

def test_init_from_seq():
    x = Chromosome('chr1',['a','c','g','t'])
    assert True

def test_slice(chr1):
    x = chr1[1:100]
    assert len(x) == 100
    
def test_invalid_slice(chr1):
    with pytest.raises(ValueError):
        chr1[0:100]

def test_get_bp(chr1):
    assert chr1[12] == 'A' 

def test_get_bp_invalid_coordinate(chr1):
    with pytest.raises(ValueError):
        chr1[0]

def test_repr(chr1):
    assert repr(chr1) == "Chromosome('AAAAAAAAAAAA...AAAAAAAAAAAAA')"

def test_equals(chr1,chr2):
    x = chr1
    y = chr2
    assert x == x
    assert x != y
