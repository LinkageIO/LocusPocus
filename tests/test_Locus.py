import pytest
import numpy as np

from itertools import chain
from locuspocus import Locus

@pytest.fixture
def simple_Locus():
    return Locus(1,100,200) 

def test_initialization(simple_Locus):
    # numeric chromosomes
    assert simple_Locus.chrom == '1'
    assert simple_Locus.start == 100
    assert simple_Locus.end == 200
    assert len(simple_Locus) == 101

def test_as_dict(simple_Locus):
    simple_Locus_dict = {
            'name': None,
            'chrom': '1',
            'start': 100,
            'end': 200,
            'feature_type': 'locus',
            'frame': None,
            'score': np.nan,
            'source': 'locuspocus',
            'strand': None,
            'window': 0
            }
    assert simple_Locus.as_dict() == simple_Locus_dict

def test_update(simple_Locus):
    new_data = {'name': 'simple_Locus'}
    new_Locus = simple_Locus.update(new_data)
    assert set(new_data.keys()) <= set(new_Locus.attr.keys()) 

    # probably we should also check the values on a key per key basis
    for k in new_data.keys():
        assert new_data[k] == new_Locus[k]

#def test_as_record(simple_Locus):
#    assert simple_Locus.as_record() == ('1', 100, 200, '<None>1:100-200', 0, 0, 0, 'LocusPocus', '<None>1:100-200')

#def test_from_record(simple_Locus):
#    assert Locus.from_record(
#        ('1', 100, 200, '<None>1:100-200', 0, 0, 0, '<None>1:100-200')) == simple_Locus

def test_setitem(simple_Locus):
    simple_Locus.attr['name'] = 'simple_Locus'
    assert simple_Locus['name'] == 'simple_Locus'

def test_getitem(simple_Locus):
    simple_Locus.attr['name'] = 'simple_Locus'
    assert simple_Locus['name'] == simple_Locus.attr['name']

#def test_update(simple_Locus):
#    assert False

def test_default_getitem(simple_Locus):
    assert simple_Locus.default_getitem('name', 'default') == 'default'

def test_start(simple_Locus):
    assert simple_Locus.start == 100

def test_end(simple_Locus):
    assert simple_Locus.end == 200

def test_coor(simple_Locus):
    assert simple_Locus.coor == (100, 200)

def test_upstream(simple_Locus):
    assert simple_Locus.upstream == 100

def test_downstream(simple_Locus):
    assert simple_Locus.downstream == 200

def test_name(simple_Locus):
    assert simple_Locus.name is None

def test_add_same_chrom(simple_Locus):
    another_Locus = Locus(1, 110, 220)
    assert simple_Locus + another_Locus == Locus(1, 100, 220)

def test_add_diff_chrom(simple_Locus):
    another_Locus = Locus(2, 110, 220)
    with pytest.raises(TypeError):
        simple_Locus + another_Locus

def test_eq(simple_Locus):
    another_Locus = Locus(1, 110, 220)
    assert simple_Locus == simple_Locus
    assert simple_Locus != another_Locus

def test_contains(simple_Locus):
    overlap_Locus = Locus(1, 110, 220)
    disjoint_Locus = Locus(1, 210, 300)
    assert simple_Locus in overlap_Locus
    assert simple_Locus not in disjoint_Locus

def test_len(simple_Locus):
    assert len(simple_Locus) == 101
    assert len(Locus(1, 100, 100)) == 1

#def test_cmp(simple_Locus):
#    assert False

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

def test_sub_same_chrom(simple_Locus):
    another_Locus = Locus(1, 110, 220)
    assert simple_Locus - another_Locus == -90

def test_sub_diff_chrom(simple_Locus):
    another_Locus = Locus(2, 110, 220)
    assert simple_Locus - another_Locus == float('inf')

def test_str(simple_Locus):
    assert str(simple_Locus) == '<None>1:100-200+0(0)'

def test_repr(simple_Locus):
    assert repr(simple_Locus) == '<None>1:100-200+0(0)'

def test_hash(simple_Locus):
    assert hash(simple_Locus) == 1035530407970800134

def test_summary(simple_Locus):
    assert simple_Locus.summary()

#def test_candidate_vs_bootstrap_length(testRefGen,testGWAS):
#    Term = next(testGWAS.iter_terms())
#    snps = Term.effective_loci(window_size=50000)
#    candidates = testRefGen.candidate_genes(snps,chain=False)
#    bootstraps = testRefGen.bootstrap_candidate_genes(snps,chain=False)
#    # Make sure we are pulling out the same number of random genes for
#    # Each locus
#    for c,b in zip(candidates,bootstraps):
#        assert len(c) == len(b)
#    assert len(set(chain(*candidates))) == len(set(chain(*bootstraps)))

#def test_generate_from_id(Zm5bFGS):
#   random_gene = Zm5bFGS.random_gene()
#   assert random_gene == Zm5bFGS[random_gene.id]


