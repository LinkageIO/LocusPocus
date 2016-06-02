import pytest

from itertools import chain
from locuspocus import Locus

@pytest.fixture
def simple_Locus():
    return Locus(1,100,200) 

def test__initialization(simple_Locus):
    # numeric chromosomes
    assert simple_Locus.chrom == '1'
    assert simple_Locus.start == 100
    assert simple_Locus.end == 200
    assert len(simple_Locus) == 100

def test_as_dict(simple_Locus):
    assert False

def test_id_getter(simple_Locus):
    assert False

def test_update(simple_Locus):
    assert False

def test_as_record(simple_Locus):
    assert False

def test_from_record(simple_Locus):
    assert False

def test_setitem(simple_Locus):
    assert False

def test_getitem(simple_Locus):
    assert False

def test_update(simple_Locus):
    assert False

def test_default_getitem(simple_Locus):
    assert False

def test_start(simple_Locus):
    assert False

def test_end(simple_Locus):
    assert False

def test_coor(simple_Locus):
    assert False

def test_upstream(simple_Locus):
    assert False

def test_downstream(simple_Locus):
    assert False

def test_name(simple_Locus):
    assert False

def test_add(simple_Locus):
    assert False

def test_eq(simple_Locus):
    assert False

def test_contains(simple_Locus):
    assert False

def test_len(simple_Locus):
    assert False

def test_cmp(simple_Locus):
    assert False

def test_lt(simple_Locus):
    assert False

def test_gt(simple_Locus):
    assert False

def test_sub(simple_Locus):
    assert False

def test_str(simple_Locus):
    assert False

def test_repr(simple_Locus):
    assert False

def test_hash(simple_Locus):
    assert False

def test_summary(simple_Locus):
    assert False

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


