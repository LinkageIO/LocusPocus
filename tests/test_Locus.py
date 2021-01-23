import pytest
import numpy as np

from itertools import chain
from locuspocus import Locus

from locuspocus.exceptions import StrandError, ChromosomeError


@pytest.fixture
def simple_Locus():
    return Locus("1", 100, 200, attrs={"foo": "bar"})


def test_initialization(simple_Locus):
    # numeric chromosomes
    assert simple_Locus.chromosome == "1"
    assert simple_Locus.start == 100
    assert simple_Locus.end == 200
    assert len(simple_Locus) == 101


def test_getitem(simple_Locus):
    assert simple_Locus["foo"] == "bar"


def test_default_getitem(simple_Locus):
    assert simple_Locus.default_getitem("name", "default") == "default"


def test_start(simple_Locus):
    assert simple_Locus.start == 100


def test_plus_stranded_start():
    l = Locus("1", 1, 100, strand="+")
    assert l.stranded_start == 1


def test_minus_stranded_start():
    l = Locus("1", 1, 100, strand="-")
    assert l.stranded_start == 100


def test_end(simple_Locus):
    assert simple_Locus.end == 200


def test_plus_stranded_end():
    l = Locus("1", 1, 100, strand="+")
    assert l.stranded_end == 100


def test_minus_stranded_end():
    l = Locus("1", 1, 100, strand="-")
    assert l.stranded_end == 1


def test_hash():
    l = Locus("1", 1, 100, strand="+")
    assert hash(l) == 530409172339088127


def test_coor(simple_Locus):
    assert simple_Locus.coor == (100, 200)


def test_upstream(simple_Locus):
    assert simple_Locus.upstream(50) == 50


def test_upstream_minus_strand(simple_Locus):
    l = Locus("1", 1, 100, strand="-")
    assert l.upstream(50) == 150


def test_downstream(simple_Locus):
    assert simple_Locus.downstream(50) == 250


def test_downstream_minus_strand(simple_Locus):
    l = Locus("1", 100, 200, strand="-")
    assert l.downstream(50) == 50


def test_center():
    l = Locus("1", 100, 200, strand="-")
    assert l.center == 150.5


def test_name(simple_Locus):
    assert simple_Locus.name is None


def test_eq(simple_Locus):
    another_Locus = Locus(1, 110, 220)
    assert simple_Locus == simple_Locus
    assert simple_Locus != another_Locus


def test_ne_diff_attrs():
    x = Locus(1, 110, 220, attrs={"foo": "bar"})
    y = Locus(1, 110, 220, attrs={"baz": "bat"})
    assert x != y


def test_empty_subloci_getitem():
    a = Locus("1", 10, 20)
    with pytest.raises(IndexError):
        a.subloci[1]


def test_empty_subloci_repr():
    a = Locus("1", 10, 20)
    b = Locus("1", 20, 30)
    x = Locus(1, 110, 220, attrs={"foo": "bar"}, subloci=[a, b])
    assert repr(x)


def test_ne_diff_subloci():
    a = Locus("1", 10, 20)
    b = Locus("1", 20, 30)
    c = Locus("1", 30, 40)
    d = Locus("1", 40, 50)

    x = Locus(1, 110, 220, attrs={"foo": "bar"}, subloci=[a, b])
    y = Locus(1, 110, 220, attrs={"foo": "bar"}, subloci=[c, d])
    assert x != y


def test_loci_lt_by_chrom():
    x = Locus("1", 1, 1)
    y = Locus("2", 1, 1)
    assert x < y


def test_loci_gt_by_chrom():
    x = Locus("1", 1, 1)
    y = Locus("2", 1, 1)
    assert y > x


def test_loci_lt_by_pos():
    x = Locus("1", 1, 100)
    y = Locus("1", 2, 100)
    assert x < y


def test_loci_gt_by_pos():
    x = Locus("1", 1, 100)
    y = Locus("1", 2, 200)
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


def test_subloci_repr(simple_Locus):
    assert repr(simple_Locus.subloci)


def test_subloci_getitem():
    x = Locus("1", 1, 2)
    y = Locus("1", 3, 4, name="sublocus")
    x.add_sublocus(y)
    assert x.subloci[0].name == "sublocus"


def test_subloci_iter():
    x = Locus("1", 1, 2)
    y = Locus("1", 3, 4, name="sublocus1")
    z = Locus("1", 3, 4, name="sublocus2")

    x.add_sublocus(y)
    x.add_sublocus(z)
    for sub in x.subloci:
        assert sub.chromosome == "1"


def test_subloci_len():
    x = Locus("1", 1, 2)
    y = Locus("1", 3, 4, name="sublocus1")
    z = Locus("1", 3, 4, name="sublocus2")

    x.add_sublocus(y)
    x.add_sublocus(z)
    assert len(x.subloci) == 2


def test_empty_attr():
    # empty attr
    x = Locus("chr1", 1, 100)
    assert x.attrs.keys() == []


def test_empty_attr_cmp():
    x = Locus("chr1", 1, 100)
    y = Locus("chr1", 2, 200)
    # both are empty
    assert x.attrs == y.attrs


def test_empty_attrs_contains():
    x = Locus("chr1", 1, 100)
    assert "foo" not in x.attrs


def test_empty_attr_values():
    # empty attr
    x = Locus("chr1", 1, 100)
    assert x.attrs.values() == []


def test_empty_getitem():
    # empty attr
    x = Locus("chr1", 1, 100)
    with pytest.raises(KeyError):
        x["test"]


def test_empty_items():
    # empty attr
    x = Locus("chr1", 1, 100)
    assert x.attrs.items() == {}


def test_attrs_repr():
    x = Locus("chr1", 1, 100, attrs={"foo": "bar"})
    assert repr(x)


def test_attrs_cmp_with_empty():
    x = Locus("chr1", 1, 100, attrs={"foo": "bar"})
    y = Locus("chr1", 2, 200)
    assert x.attrs != y.attrs


def test_empty_setitem():
    # empty attr
    x = Locus("chr1", 1, 100)
    x["test"] = "foo"


def test_attrs_keys_method():
    from locuspocus.locus import LocusAttrs

    x = LocusAttrs(attrs={"foo": "locus1", "bar": "baz"})
    assert sorted(x.keys()) == ["bar", "foo"]


def test_attrs_keys_method_empty():
    x = Locus("1", 3, 4, attrs={})
    assert len(list(x.attrs.keys())) == 0


def test_attrs_vals_method():
    from locuspocus.locus import LocusAttrs

    x = LocusAttrs(attrs={"foo": "locus1", "bar": "baz"})
    assert len(sorted(x.values())) == 2


def test_attrs_vals_method_empty():
    x = Locus("1", 3, 4, attrs={})
    assert len(list(x.attrs.values())) == 0


def test_attrs_getitem():
    from locuspocus.locus import LocusAttrs

    x = LocusAttrs(attrs={"foo": "locus1", "bar": "baz"})
    assert x["foo"] == "locus1"


def test_attrs_getitem_missing():
    from locuspocus.locus import LocusAttrs

    x = LocusAttrs(attrs={"foo": "locus1", "bar": "baz"})
    with pytest.raises(KeyError):
        x["foobar"]


def test_attrs_setitem():
    from locuspocus.locus import LocusAttrs

    x = LocusAttrs(attrs={"foo": "locus1", "bar": "baz"})
    assert x["foo"] == "locus1"
    x["foo"] = "bar"
    assert x["foo"] == "bar"


def test_setitem():
    x = Locus("1", 3, 4, attrs={"foo": "locus1", "bar": "baz"})
    assert x["foo"] == "locus1"
    x["foo"] = "bar"
    assert x["foo"] == "bar"


def test_attrs_items():
    from locuspocus.locus import LocusAttrs

    x = LocusAttrs(attrs={"foo": "locus1", "bar": "baz"})
    assert len(sorted(x.items())) == 2


def test_attrs_contains():
    x = Locus("1", 3, 4, attrs={"foo": "locus1", "bar": "baz"})
    assert "foo" in x.attrs


def test_le_equals():
    x = Locus("1", 3, 4, attrs={"foo": "locus1", "bar": "baz"})
    y = Locus("1", 3, 4)
    assert x <= y


def test_le_less():
    x = Locus("1", 3, 4, attrs={"foo": "locus1", "bar": "baz"})
    y = Locus("1", 30, 40)
    assert x <= y


def test_ge_equals():
    x = Locus("1", 3, 4, attrs={"foo": "locus1", "bar": "baz"})
    y = Locus("1", 3, 4)
    assert x >= y


def test_ge_greater():
    x = Locus("1", 30, 40, attrs={"foo": "locus1", "bar": "baz"})
    y = Locus("1", 3, 4)
    assert x >= y


def test_stranded_start_invalid():
    # Strand cannot be '='
    x = Locus("1", 3, 4, strand="=")
    with pytest.raises(StrandError):
        x.stranded_start


def test_stranded_stop_invalid():
    # Strand cannot be '='
    x = Locus("1", 3, 4, strand="=")
    with pytest.raises(StrandError):
        x.stranded_end


def test_as_record():
    x = Locus("1", 3, 4, strand="+")
    # This doesn't compare the dictionaries of each ...
    assert x.as_record()[0] == (
        "1",
        3,
        4,
        "locuspocus",
        "locus",
        "+",
        None,
        None,
        2039807104618252476,
    )


def test_center_distance():
    x = Locus("1", 1, 100, strand="+")
    # This needs to be 201 since x starts at 1
    y = Locus("1", 201, 300, strand="=")
    assert x.center_distance(y) == 200


def test_center_distance_different_chroms():
    x = Locus("1", 1, 100, strand="+")
    # This needs to be 201 since x starts at 1
    y = Locus("2", 201, 300, strand="+")
    assert x.center_distance(y) == np.inf


def test_str():
    x = Locus("1", 1, 100, strand="+")
    assert str(x) == repr(x)


def test_combine():
    x = Locus("1", 1, 2)
    y = Locus("1", 3, 4)
    z = x.combine(y)

    assert z.start == 1
    assert z.end == 4
    assert x in z.subloci
    assert y in z.subloci


def test_combine_chromosome_mismatch():
    x = Locus("1", 1, 2)
    y = Locus("2", 3, 4)
    with pytest.raises(ChromosomeError):
        x.combine(y)


def test_distance():
    x = Locus("1", 1, 100)
    y = Locus("1", 150, 250)
    assert x.distance(y) == 49


def test_distance_diff_chroms():
    x = Locus("1", 1, 100)
    y = Locus("2", 150, 250)
    assert x.distance(y) == np.inf
