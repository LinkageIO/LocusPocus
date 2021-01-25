

"""
    Unit tests for Ontology
"""

import minus80 as m80
import locuspocus as lp

def test_init(testOnt):
    try:
        testOnt
        return True
    except NameError:
        return False

def test_len(testOnt):
    assert len(testOnt) == 10


def test_iter(testOnt):
    for term in testOnt:
        assert isinstance(term, lp.Term)

def test_access_loci(testOnt):
    assert isinstance(testOnt.loci, lp.Loci)

def test_add_term():
    try:
        x = lp.Ontology("empty")
        assert len(x) == 0
        x.add_term(lp.Term("test",loci=[lp.Locus(1,1,1)]))
        assert len(x) == 1
    finally:
        m80.delete("Ontology","empty")

def test_num_terms(testOnt):
    assert testOnt.num_terms() == len(testOnt)

def test_num_terms_with_min(testOnt):
    # should match 6,7,8,9,10
    assert testOnt.num_terms(min_term_size=6) == 5

def test_num_terms_with_max(testOnt):
    # should match 1,2,3,4,5
    assert testOnt.num_terms(max_term_size=5) == 5

def test_get_item(testOnt):
    assert isinstance(testOnt['term_1'], lp.Term)

def test_get_item_by_TID(testOnt):
    # TIDs start at 1
    assert isinstance(testOnt[1], lp.Term)

def test_terms_containing(testOnt):
    # Locus 1,1,1 should be in all terms
    assert len(testOnt.terms_containing([lp.Locus(1,1,1)])) == len(testOnt) 

def test_terms_function(testOnt):
    for term in testOnt.terms():
        assert isinstance(term, lp.Term)

def test_terms_with_min(testOnt):
    # should match 6,7,8,9,10
    assert len(list(testOnt.terms(min_term_size=6))) == 5

def test_terms_with_max(testOnt):
    # should match 1,2,3,4,5
    assert len(list(testOnt.terms(max_term_size=5))) == 5

def test_summary(testOnt):
    assert isinstance(testOnt.summary(),str)

def test_rand(testOnt):
    assert isinstance(testOnt.rand(), lp.Term)
