import os
import pytest

from locuspocus import Locus, Loci
from locuspocus.exceptions import MissingLocusError, StrandError

import minus80 as m80

"""
    Unit tests for Loci
"""

NUM_GENES = 39656


def test_init(testRefGen):
    try:
        testRefGen
        return True
    except NameError:
        return False


def test_len(testRefGen):
    "quick check to make sure that the refgen reports the correct number of features"
    assert len(testRefGen) == NUM_GENES


def test_get_locus_by_LID(testRefGen):
    "make sure that fetching a locus by its LID yields the same locus"
    rand_locus = testRefGen.rand()
    assert rand_locus == testRefGen._get_locus_by_LID(rand_locus._LID)


def test_get_locus_by_LID_missing(testRefGen):
    "make sure that fetching a locus by its LID yields the same locus"
    with pytest.raises(MissingLocusError):
        testRefGen._get_locus_by_LID(-1)


def test_get_LID_missing(testRefGen):
    "Make sure that fetching a locus by LID returns the same locus"
    with pytest.raises(MissingLocusError):
        assert testRefGen._get_LID(Locus("na", 1, 1))


def test_get_LID_from_name_missing(testRefGen):
    "Make sure that fetching a locus by LID returns the same locus"
    with pytest.raises(MissingLocusError):
        assert testRefGen._get_LID("DoesNotExist")


def test_add_locus():
    "add a locus to an empty refloci db and then retrieve it"
    if m80.exists("Loci", "empty"):
        m80.delete("Loci", "empty")
    empty = Loci("empty")
    assert len(empty) == 0
    empty.add_locus(Locus("1", 1, 1, feature_type="gene", attrs={"foo": "bar"}))
    assert len(empty) == 1
    m80.delete("Loci", "empty")


def test_add_locus_with_attrs():
    "add a locus to an empty refloci db and then retrieve it"
    if m80.exists("Loci", "empty"):
        m80.delete("Loci", "empty")
    empty = Loci("empty")
    assert len(empty) == 0
    LID = empty.add_locus(Locus("1", 1, 1, feature_type="gene", attrs={"foo": "bar"}))
    assert len(empty) == 1
    l = empty._get_locus_by_LID(LID)
    assert l["foo"] == "bar"
    m80.delete("Loci", "empty")


def test_nuke_tables():
    "add a locus to an empty refloci db and then retrieve it"
    if m80.exists("Loci", "empty"):
        m80.delete("Loci", "empty")
    empty = Loci("empty")
    assert len(empty) == 0
    x = Locus("1", 1, 1, feature_type="gene", attrs={"foo": "bar"})
    y = Locus("1", 2, 2, feature_type="exon", attrs={"baz": "bat"})
    x.add_sublocus(y)
    empty.add_locus(x)
    assert len(empty) == 1
    empty._nuke_tables()
    assert len(empty) == 0
    m80.delete("Loci", "empty")


def test_add_locus_with_subloci():
    "add a locus to an empty refloci db and then retrieve it"
    if m80.exists("Loci", "empty"):
        m80.delete("Loci", "empty")
    empty = Loci("empty")
    assert len(empty) == 0
    x = Locus("1", 1, 1, feature_type="gene", attrs={"foo": "bar"})
    y = Locus("1", 2, 2, feature_type="exon", attrs={"baz": "bat"})
    x.add_sublocus(y)
    LID = empty.add_locus(x)
    assert len(empty) == 1
    l = empty._get_locus_by_LID(LID)
    assert l["foo"] == "bar"
    assert len(l.subloci) == 1
    m80.delete("Loci", "empty")


def test_import_gff(testRefGen):
    "test importing loci from a GFF file"
    # as the testRefGen fixture is built from a GFF
    # this will only pass if it is built
    assert testRefGen


# def test_contains_true(testRefGen):
#    'get a random locus and then test it is in the Loci object'
#    assert testRefGen.rand() in testRefGen


def test_contains_false(testRefGen):
    assert ("NO" in testRefGen) is False


def test_get_item(testRefGen):
    """
    Get a random locus and then
    retrieve that locus again
    by its id
    """
    random_locus = testRefGen.rand()
    assert random_locus == testRefGen[random_locus.name]


def test_iter(testRefGen):
    "test that the iter interface works"
    i = 0
    for locus in testRefGen:
        i += 1
    assert i == NUM_GENES


def test_rand(testRefGen):
    "test instance type"
    assert isinstance(testRefGen.rand(), Locus)


def test_rand_length(testRefGen):
    "test the length of rand with n specified"
    assert len(testRefGen.rand(n=100)) == 100


def test_rand_distinct(testRefGen):
    assert len(testRefGen.rand(2000, distinct=True)) == 2000


def test_rand_too_many(testRefGen):
    try:
        testRefGen.rand(100000)
    except ValueError:
        assert True


def test_rand_no_autopop(testRefGen):
    assert len(testRefGen.rand(1, autopop=False)) == 1


# The first 4 genes on chromosome 9
# 1       ensembl gene    4854    9652    .       -       .       ID=GRMZM2G059865;Name=GRMZM2G059865;biotype=protein_coding
# 1       ensembl gene    9882    10387   .       -       .       ID=GRMZM5G888250;Name=GRMZM5G888250;biotype=protein_coding
# 1       ensembl gene    109519  111769  .       -       .       ID=GRMZM2G093344;Name=GRMZM2G093344;biotype=protein_coding
# 1       ensembl gene    136307  138929  .       +       .       ID=GRMZM2G093399;Name=GRMZM2G093399;biotype=protein_coding


def test_within(testRefGen):
    "simple within to get chromosomal segment"
    assert len(list(testRefGen.within(Locus("1", 1, 139000), partial=False))) == 4


def test_within_bad_strand(testRefGen):
    "simple within to get chromosomal segment"
    with pytest.raises(StrandError):
        assert (
            len(
                list(
                    testRefGen.within(Locus("1", 1, 139000, strand="="), partial=False)
                )
            )
            == 4
        )


def test_within_yields_nothing(testRefGen):
    l = Locus("0", 1, 1, strand="+")
    assert len(list(testRefGen.within(l, partial=False))) == 0


def test_within_partial_false(testRefGen):
    "put the locus boundaries within gene [1] and [4] and exclude them with partial"
    assert len(list(testRefGen.within(Locus("1", 6000, 137000), partial=False))) == 2


def test_within_partial_true(testRefGen):
    "put the locus boundaries within gene [1] and [4] and exclude them with partial"
    assert len(list(testRefGen.within(Locus("1", 6000, 137000), partial=True))) == 4


def test_within_same_strand(testRefGen):
    "test fetching loci only on the same strand"
    assert (
        len(
            list(
                testRefGen.within(
                    Locus("1", 1, 139000, strand="+"), partial=True, same_strand=True
                )
            )
        )
        == 1
    )


def test_within_same_strand_and_ignore_strand(testRefGen):
    "test fetching loci only on the same strand"
    with pytest.raises(ValueError):
        list(
            testRefGen.within(
                Locus("1", 1, 139000, strand="+"), ignore_strand=True, same_strand=True
            )
        )


def test_within_same_strand_minus(testRefGen):
    "test fetching loci only on the same strand"
    assert (
        len(
            list(
                testRefGen.within(
                    Locus("1", 1, 139000, strand="-"), partial=True, same_strand=True
                )
            )
        )
        == 3
    )


def test_within_strand_order(testRefGen):
    # should return the locus at the beginning of the chromosome
    loci = list(testRefGen.within(Locus("1", 1, 139000, strand="+")))
    assert loci[0].start == 4854


def test_within_strand_order_minus(testRefGen):
    # should return fourth gene first
    loci = list(testRefGen.within(Locus("1", 1, 139000, strand="-")))
    assert loci[0].start == 136307


def test_within_error_on_both_same_strand_and_ignore_strand(testRefGen):
    try:
        testRefGen.within(
            Locus("1", 1, 139000, strand="-"), ignore_strand=True, same_strand=True
        )
    except ValueError:
        assert True


def test_upstream_plus_strand(testRefGen):
    # Below is GRMZM2G093399, but on the minus strand
    x = Locus("1", 136307, 138929)
    l = [x.name for x in testRefGen.upstream_loci(x, n=3)]
    assert l[0] == "GRMZM2G093344"
    assert l[1] == "GRMZM5G888250"
    assert l[2] == "GRMZM2G059865"


def test_upstream_minus_strand(testRefGen):
    x = testRefGen["GRMZM5G888250"]
    l = [x.name for x in testRefGen.upstream_loci(x, n=2)]
    assert l[0] == "GRMZM2G093344"
    assert l[1] == "GRMZM2G093399"


def test_upstream_accepts_loci(testRefGen):
    loci = [testRefGen["GRMZM2G093399"], testRefGen["GRMZM2G093399"]]
    l1, l2 = map(list, testRefGen.upstream_loci(loci, n=2))
    assert len(l1) == 2
    assert len(l2) == 2


def test_upstream_same_strand(testRefGen):
    x = testRefGen["GRMZM5G888250"]
    for x in testRefGen.upstream_loci(x, n=5, same_strand=True):
        assert x.strand == "-"


def test_upstream_limit_n(testRefGen):
    g = testRefGen.upstream_loci(testRefGen["GRMZM2G093399"], n=2)
    assert len(list(g)) == 2


def test_upstream_filter_same_strand(testRefGen):
    g = testRefGen.upstream_loci(testRefGen["GRMZM2G093399"], n=3, same_strand=True)
    assert len(list(g)) == 0


def test_downstream(testRefGen):
    l = [x.name for x in testRefGen.downstream_loci(Locus("1", 4854, 9652), n=3)]
    assert l[0] == "GRMZM5G888250"
    assert l[1] == "GRMZM2G093344"
    assert l[2] == "GRMZM2G093399"


def test_downstream_accepts_loci(testRefGen):
    x = Locus("1", 4854, 9652)
    loci = [x, x]
    l1, l2 = map(list, testRefGen.downstream_loci(loci, n=2))
    assert len(l1) == 2
    assert len(l2) == 2


def test_downstream_limit_n(testRefGen):
    l = [x.name for x in testRefGen.downstream_loci(Locus("1", 4854, 9652), n=1)]
    assert l[0] == "GRMZM5G888250"
    assert len(l) == 1


def test_downstream_same_strand_limit_n(testRefGen):
    l = [
        x.name
        for x in testRefGen.downstream_loci(
            Locus("1", 4854, 9652), n=1, same_strand=True
        )
    ]
    assert len(l) == 1
    # This skips GRMZM5G888250 and GRMZM2G093344 since they are on the - strand
    assert l[0] == "GRMZM2G093399"


def test_flank_loci_limited_n(testRefGen):
    x = Locus("1", 10500, 10500)
    up, down = testRefGen.flanking_loci(x, n=2)
    assert len(list(up)) == 2
    assert len(list(down)) == 2


def test_encompassing_loci(testRefGen):
    x = Locus("1", 10000, 10000)
    loci = list(testRefGen.encompassing_loci(x))
    assert loci[0].name == "GRMZM5G888250"


def test_full_import_gff():
    if m80.exists("Loci", "ZmSmall"):
        m80.delete("Loci", "ZmSmall")
    gff = os.path.expanduser(os.path.join("raw", "maize_small.gff"))
    x = Loci("ZmSmall")
    x.import_gff(gff)
    m80.delete("Loci", "ZmSmall")


def test_import_gff_gzipped():
    if m80.exists("Loci", "ZmSmall"):
        m80.delete("Loci", "ZmSmall")
    gff = os.path.expanduser(os.path.join("raw", "maize_small.gff.gz"))
    x = Loci("ZmSmall")
    x.import_gff(gff)
    m80.delete("Loci", "ZmSmall")
