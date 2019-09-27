
import os
import pytest

from locuspocus import Locus, FrozenLoci
from locuspocus.exceptions import *

import minus80 as m80

'''
    Unit tests for FrozenLoci
'''

NUM_GENES = 39656
NUM_CDS = 292257

def test_init(testFrozenLoci):
    assert testFrozenLoci

def test_primary_LIDS(testFrozenLoci):
    'make sure the internal _primary_LIDS method works'
    with testFrozenLoci.filter_feature_type('gene'):
        assert len(testFrozenLoci._primary_LIDS()) == NUM_GENES

def test_len(testFrozenLoci):
    'quick check to make sure that the refgen reports the correct number of features'
    with testFrozenLoci.filter_feature_type('gene'):
        assert len(testFrozenLoci) == NUM_GENES

def test_get_locus_by_LID(testFrozenLoci):
    'make sure that fetching a locus by its LID yields the same locus'
    rand_locus = testFrozenLoci.rand()
    assert rand_locus == testFrozenLoci._get_locus_by_LID(rand_locus._LID)

def test_get_locus_by_LID_missing(testFrozenLoci):
    'make sure that fetching a locus by its LID yields the same locus'
    with pytest.raises(MissingLocusError):
        testFrozenLoci._get_locus_by_LID(-1)

#def test_get_LID(testFrozenLoci):
#    'Make sure that fetching a locus by LID returns the same locus'
#    x = testFrozenLoci.rand()
#    assert testFrozenLoci._get_LID(x) == x._LID

def test_get_LID_missing(testFrozenLoci):
    'Make sure that fetching a locus by LID returns the same locus'
    with pytest.raises(MissingLocusError):
        assert testFrozenLoci._get_LID(Locus('na',1,1))

def test_get_LID_from_name_missing(testFrozenLoci):
    'Make sure that fetching a locus by LID returns the same locus'
    with pytest.raises(MissingLocusError):
        assert testFrozenLoci._get_LID('DoesNotExist')

def test_add_locus():
    'add a locus to an empty refloci db and then retrieve it'
    if m80.Tools.available('FrozenLoci','empty'):
        m80.Tools.delete('FrozenLoci','empty',force=True)
    empty = FrozenLoci('empty')
    assert len(empty) == 0
    empty.add_locus(Locus('1',1,1,feature_type='gene',attrs={'foo':'bar'}))
    assert len(empty) == 1
    m80.Tools.delete('FrozenLoci','empty',force=True)

def test_add_locus_with_attrs():
    'add a locus to an empty refloci db and then retrieve it'
    if m80.Tools.available('FrozenLoci','empty'):
        m80.Tools.delete('FrozenLoci','empty',force=True)
    empty = FrozenLoci('empty')
    assert len(empty) == 0
    LID = empty.add_locus(Locus('1',1,1,feature_type='gene',attrs={'foo':'bar'}))
    assert len(empty) == 1
    l = empty._get_locus_by_LID(LID)
    assert l['foo'] == 'bar'
    m80.Tools.delete('FrozenLoci','empty',force=True)

def test_nuke_tables():
    'add a locus to an empty refloci db and then retrieve it'
    if m80.Tools.available('FrozenLoci','empty'):
        m80.Tools.delete('FrozenLoci','empty',force=True)
    empty = FrozenLoci('empty')
    assert len(empty) == 0
    x = Locus('1',1,1,feature_type='gene',attrs={'foo':'bar'})
    y = Locus('1',2,2,feature_type='exon',attrs={'baz':'bat'})
    x.add_sublocus(y)
    LID = empty.add_locus(x)
    assert len(empty) == 1
    empty._nuke_tables()
    assert len(empty) == 0
    m80.Tools.delete('FrozenLoci','empty',force=True)


def test_add_locus_with_subloci():
    'add a locus to an empty refloci db and then retrieve it'
    if m80.Tools.available('FrozenLoci','empty'):
        m80.Tools.delete('FrozenLoci','empty',force=True)
    empty = FrozenLoci('empty')
    assert len(empty) == 0
    x = Locus('1',1,1,feature_type='gene',attrs={'foo':'bar'})
    y = Locus('1',2,2,feature_type='exon',attrs={'baz':'bat'})
    x.add_sublocus(y)
    LID = empty.add_locus(x)
    assert len(empty) == 1
    l = empty._get_locus_by_LID(LID)
    assert l['foo'] == 'bar'
    assert len(l.subloci) == 1
    m80.Tools.delete('FrozenLoci','empty',force=True)


def test_import_gff(testFrozenLoci):
    'test importing loci from a GFF file'
    # as the testFrozenLoci fixture is built from a GFF
    # this will only pass if it is built
    assert testFrozenLoci

#def test_contains_true(testFrozenLoci):
#    'get a random locus and then test it is in the FrozenLoci object'
#    assert testFrozenLoci.rand() in testRefGen

def test_contains_false(testFrozenLoci):
    assert ('NO' in testFrozenLoci) is False

def test_get_item(testFrozenLoci):
    '''
        Get a random locus and then
        retrieve that locus again
        by its id
    '''
    random_locus = testFrozenLoci.rand()
    assert random_locus == testFrozenLoci[random_locus.name]

def test_iter(testFrozenLoci):
    'test that the iter interface works' 
    i = 0
    for locus in testFrozenLoci:
        i += 1
    assert i == NUM_GENES

def test_filter_feature_type_context_manager(testFrozenLoci):
    # Check to see we are in gene space
    prev_len = len(testFrozenLoci)
    # temporarily conver to CDS
    with testFrozenLoci.filter_feature_type('CDS'):
        assert len(testFrozenLoci) == NUM_CDS
    # assert we didn't need to do anything to switch back
    assert len(testFrozenLoci) == prev_len

def test_set_primary_feature_type(testFrozenLoci):
    testFrozenLoci.set_primary_feature_type('gene')
    assert len(testFrozenLoci) == NUM_GENES
    testFrozenLoci.set_primary_feature_type('CDS')
    assert len(testFrozenLoci) == NUM_CDS
    # switch back to gene
    testFrozenLoci.set_primary_feature_type('gene')
    assert len(testFrozenLoci) == NUM_GENES

def test_rand(testFrozenLoci):
    'test instance type'
    assert isinstance(testFrozenLoci.rand(),Locus)

def test_rand_length(testFrozenLoci):
    'test the length of rand with n specified'
    assert len(testFrozenLoci.rand(n=100)) == 100 

def test_rand_distinct(testFrozenLoci):
    assert len(testFrozenLoci.rand(2000,distinct=True)) == 2000

def test_rand_too_many(testFrozenLoci):
    try:
        testFrozenLoci.rand(100000)
    except ValueError as e:
        assert True

def test_rand_no_autopop(testFrozenLoci):
    assert len(testFrozenLoci.rand(1,autopop=False)) == 1

# The first 4 genes on chromosome 9
# 1       ensembl gene    4854    9652    .       -       .       ID=GRMZM2G059865;Name=GRMZM2G059865;biotype=protein_coding
# 1       ensembl gene    9882    10387   .       -       .       ID=GRMZM5G888250;Name=GRMZM5G888250;biotype=protein_coding
# 1       ensembl gene    109519  111769  .       -       .       ID=GRMZM2G093344;Name=GRMZM2G093344;biotype=protein_coding
# 1       ensembl gene    136307  138929  .       +       .       ID=GRMZM2G093399;Name=GRMZM2G093399;biotype=protein_coding

def test_within(testFrozenLoci):
    'simple within to get chromosomal segment'
    assert len(list(testFrozenLoci.within(Locus('1',1,139000),partial=False))) == 4

def test_within_bad_strand(testFrozenLoci):
    'simple within to get chromosomal segment'
    with pytest.raises(StrandError):
        assert len(list(testFrozenLoci.within(Locus('1',1,139000,strand='='),partial=False))) == 4

def test_within_yields_nothing(testFrozenLoci):
    l = Locus('0',1,1,strand='+')
    assert len(list(testFrozenLoci.within(l,partial=False))) == 0

def test_within_partial_false(testFrozenLoci):
    'put the locus boundaries within gene [1] and [4] and exclude them with partial'
    assert len(list(testFrozenLoci.within(Locus('1',6000,137000),partial=False))) == 2

def test_within_partial_true(testFrozenLoci):
    'put the locus boundaries within gene [1] and [4] and exclude them with partial'
    assert len(list(testFrozenLoci.within(Locus('1',6000,137000),partial=True))) == 4

def test_within_same_strand(testFrozenLoci):
    'test fetching loci only on the same strand'
    assert len(list(testFrozenLoci.within(Locus('1',1,139000,strand='+'),partial=True,same_strand=True))) == 1

def test_within_same_strand_and_ignore_strand(testFrozenLoci):
    'test fetching loci only on the same strand'
    with pytest.raises(ValueError):
        list(testFrozenLoci.within(Locus('1',1,139000,strand='+'),ignore_strand=True,same_strand=True))

def test_within_same_strand_minus(testFrozenLoci):
    'test fetching loci only on the same strand'
    assert len(list(testFrozenLoci.within(Locus('1',1,139000,strand='-'),partial=True,same_strand=True))) == 3

def test_within_strand_order(testFrozenLoci):
    # should return the locus at the beginning of the chromosome 
    loci = list(testFrozenLoci.within(Locus('1',1,139000,strand='+')))
    assert loci[0].start == 4854

def test_within_strand_order_minus(testFrozenLoci):
    # should return fourth gene first 
    loci = list(testFrozenLoci.within(Locus('1',1,139000,strand='-')))
    assert loci[0].start == 136307

def test_within_error_on_both_same_strand_and_ignore_strand(testFrozenLoci):
    try:
        testFrozenLoci.within(Locus('1',1,139000,strand='-'),ignore_strand=True,same_strand=True)
    except ValueError as e:
        assert True

def test_upstream_plus_strand(testFrozenLoci):
    # Below is GRMZM2G093399, but on the minus strand
    x = Locus('1',136307,138929)
    l = [x.name for x in testFrozenLoci.upstream_loci(x,n=3)]
    assert l[0] == 'GRMZM2G093344'
    assert l[1] == 'GRMZM5G888250'
    assert l[2] == 'GRMZM2G059865'

def test_upstream_minus_strand(testFrozenLoci):
    x = testFrozenLoci['GRMZM5G888250']
    l = [x.name for x in testFrozenLoci.upstream_loci(x,n=2)]
    assert l[0] == 'GRMZM2G093344'
    assert l[1] == 'GRMZM2G093399'

def test_upstream_accepts_loci(testFrozenLoci):
    loci = [testFrozenLoci['GRMZM2G093399'],testRefGen['GRMZM2G093399']]
    l1,l2 = map(list,testFrozenLoci.upstream_loci(loci,n=2))
    assert len(l1) == 2
    assert len(l2) == 2

def test_upstream_same_strand(testFrozenLoci):
    x = testFrozenLoci['GRMZM5G888250']
    for x in testFrozenLoci.upstream_loci(x,n=5,same_strand=True):
        assert x.strand == '-'

def test_upstream_limit_n(testFrozenLoci):
    g = testFrozenLoci.upstream_loci(testRefGen['GRMZM2G093399'],n=2)
    assert len(list(g)) == 2

def test_upstream_filter_same_strand(testFrozenLoci):
    g = testFrozenLoci.upstream_loci(testRefGen['GRMZM2G093399'],n=3,same_strand=True)
    assert len(list(g)) == 0

def test_downstream(testFrozenLoci):
    l = [x.name for x in testFrozenLoci.downstream_loci(Locus('1',4854,9652),n=3)]
    assert l[0] == 'GRMZM5G888250'
    assert l[1] == 'GRMZM2G093344'
    assert l[2] == 'GRMZM2G093399'

def test_downstream_accepts_loci(testFrozenLoci):
    x = Locus('1',4854,9652)
    loci = [x,x]
    l1,l2 = map(list,testFrozenLoci.downstream_loci(loci,n=2))
    assert len(l1) == 2
    assert len(l2) == 2

def test_downstream_limit_n(testFrozenLoci):
    l = [x.name for x in testFrozenLoci.downstream_loci(Locus('1',4854,9652),n=1)]
    assert l[0] == 'GRMZM5G888250'
    assert len(l) == 1

def test_downstream_same_strand_limit_n(testFrozenLoci):
    l = [x.name for x in testFrozenLoci.downstream_loci(Locus('1',4854,9652),n=1,same_strand=True)]
    assert len(l) == 1
    # This skips GRMZM5G888250 and GRMZM2G093344 since they are on the - strand
    assert l[0] == 'GRMZM2G093399'

def test_flank_loci_limited_n(testFrozenLoci):
    x = Locus('1',10500,10500)
    up,down = testFrozenLoci.flanking_loci(x,n=2)
    assert len(list(up)) == 2
    assert len(list(down)) == 2

def test_encompassing_loci(testFrozenLoci):
    x = Locus('1',10000,10000)
    loci = list(testFrozenLoci.encompassing_loci(x))
    assert loci[0].name == 'GRMZM5G888250'

def test_import_gff():
    if m80.Tools.available('FrozenLoci','ZmSmall'):
        m80.Tools.delete('FrozenLoci','ZmSmall',force=True)
    gff = os.path.expanduser(
        os.path.join(
            'raw', 
            'maize_small.gff'
        )
    )
    x = FrozenLoci('ZmSmall')
    x.import_gff(gff)
    m80.Tools.delete('FrozenLoci','ZmSmall',force=True)

def test_import_gff_gzipped():
    if m80.Tools.available('FrozenLoci','ZmSmall'):
        m80.Tools.delete('FrozenLoci','ZmSmall',force=True)
    gff = os.path.expanduser(
        os.path.join(
            'raw', 
            'maize_small.gff.gz'
        )
    )
    x = FrozenLoci('ZmSmall')
    x.import_gff(gff)
    m80.Tools.delete('FrozenLoci','ZmSmall',force=True)

  
