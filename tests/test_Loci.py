import pytest

from locuspocus import Locus

'''
    Unit tests
'''

def test_from_ids(testRefGen):
    assert False

def test_get_item(testRefGen):
    assert False

def test_get_items_from_list(testRefGen):
    assert False

def test_lowercase_get_item(testRefGen):
    assert False

def test_genes_within(testRefGen):
    assert False

def test_locus_not_in_upstream_downstream(testRefGen):
    '''
        Upstream and downstream should not include the gene of interest.
    '''
    assert False

def test_upstream_downstream_genes(testRefGen):
    '''
        Take downstream of genes, then upstream genes of the 
        last gene in downstream. Tests that returns the same interval 
    '''
    assert False

def test_flanking_genes(testRefGen):
    assert False

def test_flanking_genes_includes_within_genes_for_SNPS(testRefGen):
    assert False

def test_candidate_genes_from_SNP(testRefGen):
    assert False

def test_candidate_genes_from_gene_includes_gene(testRefGen):
    assert False

def test_non_chained_candidates(testRefGen):
    assert False

def test_flank_limit_for_candidate_genes(testRefGen):
    assert False

def test_flank_limit_for_candidate_genes_from_SNP(testRefGen):
    assert False

def test_bootstrap_candidate_length_equal_from_SNP(testRefGen):
    assert False

def test_bootstrap_candidate_length_equal_from_gene(testRefGen):
    assert False

def test_refgen_length(testRefGen):
    assert False

def test_filtered_refgen(testRefGen):
    assert False


