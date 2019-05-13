import pytest

from locuspocus import Locus,RefLoci
import minus80 as m80

'''
    Unit tests
'''


def test_init(testRefGen):
    assert testRefGen

def test_primary_LIDS(testRefGen):
    'make sure the internal _primary_LIDS method works'
    with testRefGen.filter_feature_type('gene'):
        assert len(testRefGen._primary_LIDS()) == 39656

def test_len(testRefGen):
    'quick check to make sure that the refgen reports the correct number of features'
    with testRefGen.filter_feature_type('gene'):
        assert len(testRefGen) == 39656

def test_get_locus_by_LID(testRefGen):
    'make sure that fetching a locus by its LID yields the same locus'
    rand_locus = testRefGen.rand()
    assert rand_locus == testRefGen._get_locus_by_LID(rand_locus._LID)

def test_get_LID(testRefGen):
    'Make sure that fetching a locus by LID returns the same locus'
    x = testRefGen.rand()
    assert testRefGen._get_LID(x) == x._LID

def test_add_locus():
    'add a locus to an empty refloci db and then retrieve it'
    if m80.Tools.available('RefLoci','empty'):
        m80.Tools.delete('RefLoci','empty',force=True)
    empty = RefLoci('empty')
    assert len(empty) == 0
    empty.add_locus(Locus('1',1,1,feature_type='gene'))
    assert len(empty) == 1
    m80.Tools.delete('RefLoci','empty',force=True)

def test_import_gff(testRefGen):
    'test importing loci from a GFF file'
    # as the testRefGen fixture is built from a GFF
    # this will only pass if it is built
    assert testRefGen

def test_contains(testRefGen):
    'get a random locus and then test it is in the RefLoci object'
    assert testRefGen.rand() in testRefGen

def test_get_item(testRefGen):
    '''
        Get a random locus and then
        retrieve that locus again
        by its id
    '''
    random_locus = testRefGen.rand()
    assert random_locus == testRefGen[random_locus.name]

def test_iter(testRefGen):
    'test that the iter interface works' 
    i = 0
    for locus in testRefGen:
        i += 1
    assert i == 39656

def test_filter_feature_type_context_manager(testRefGen):
    # Check to see we are in gene space
    prev_len = len(testRefGen)
    # temporarily conver to CDS
    with testRefGen.filter_feature_type('CDS'):
        assert len(testRefGen) == 292257
    # assert we didn't need to do anything to switch back
    assert len(testRefGen) == prev_len

def test_set_primary_feature_type(testRefGen):
    testRefGen.set_primary_feature_type('gene')
    assert len(testRefGen) == 39656
    testRefGen.set_primary_feature_type('CDS')
    assert len(testRefGen) == 292257
    # switch back to gene
    testRefGen.set_primary_feature_type('gene')
    assert len(testRefGen) == 39656

def test_rand(testRefGen):
    'test instance type'
    assert isinstance(testRefGen.rand(),Locus)

def test_rand_length(testRefGen):
    'test the length of rand with n specified'
    assert len(testRefGen.rand(n=100)) == 100 

def test_rand_distinct(testRefGen):
    assert len(testRefGen.rand(2000,distinct=True)) == 2000

def test_rand_too_many(testRefGen):
    try:
        testRefGen.rand(100000)
    except ValueError as e:
        assert True

def test_rand_no_autopop(testRefGen):
    assert len(testRefGen.rand(1,autopop=False)) == 1

# The first 4 genes on chromosome 9
# 1       ensembl gene    4854    9652    .       -       .       ID=GRMZM2G059865;Name=GRMZM2G059865;biotype=protein_coding
# 1       ensembl gene    9882    10387   .       -       .       ID=GRMZM5G888250;Name=GRMZM5G888250;biotype=protein_coding
# 1       ensembl gene    109519  111769  .       -       .       ID=GRMZM2G093344;Name=GRMZM2G093344;biotype=protein_coding
# 1       ensembl gene    136307  138929  .       +       .       ID=GRMZM2G093399;Name=GRMZM2G093399;biotype=protein_coding


def test_within(testRefGen):
    'simple within to get chromosomal segment'
    assert len(list(testRefGen.within(Locus('1',1,139000),partial=False))) == 4

def test_within_partial_false(testRefGen):
    'put the locus boundaries within gene [1] and [4] and exclude them with partial'
    assert len(list(testRefGen.within(Locus('1',6000,137000),partial=False))) == 2

def test_within_partial_true(testRefGen):
    'put the locus boundaries within gene [1] and [4] and exclude them with partial'
    assert len(list(testRefGen.within(Locus('1',6000,137000),partial=True))) == 4

def test_within_same_strand(testRefGen):
    'test fetching loci only on the same strand'
    assert len(list(testRefGen.within(Locus('1',1,139000,strand='+'),partial=True,same_strand=True))) == 1

def test_within_same_strand_minus(testRefGen):
    'test fetching loci only on the same strand'
    assert len(list(testRefGen.within(Locus('1',1,139000,strand='-'),partial=True,same_strand=True))) == 3

def test_within_strand_order(testRefGen):
    # should return the locus at the beginning of the chromosome 
    loci = list(testRefGen.within(Locus('1',1,139000,strand='+')))
    assert loci[0].start == 4854

def test_within_strand_order_minus(testRefGen):
    # should return fourth gene first 
    loci = list(testRefGen.within(Locus('1',1,139000,strand='-')))
    assert loci[0].start == 136307

def test_within_strand_order_minus(testRefGen):
    # should return first gene since we ignore the strand 
    loci = list(testRefGen.within(Locus('1',1,139000,strand='-'),ignore_strand=True))
    assert loci[0].start == 4854

def test_within_error_on_both_same_strand_and_ignore_strand(testRefGen):
    try:
        testRefGen.within(Locus('1',1,139000,strand='-'),ignore_strand=True,same_strand=True)
    except ValueError as e:
        assert True

def test_upstream(testRefGen):
    l = [x.name for x in testRefGen.upstream_loci(testRefGen['GRMZM2G093399'],n=3)]
    assert l[0] == 'GRMZM2G093344'
    assert l[1] == 'GRMZM5G888250'
    assert l[2] == 'GRMZM2G059865'
