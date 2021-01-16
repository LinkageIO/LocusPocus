"""
    Tests
"""
import pytest
import locuspocus as lp
import minus80 as m80


def test_init(smpl_fasta):
    assert len(smpl_fasta["chr1"]) == 500000
    assert len(smpl_fasta["chr2"]) == 500000
    assert len(smpl_fasta["chr3"]) == 500000
    assert len(smpl_fasta["chr4"]) == 500000


def test_m80_init(m80_Fasta):
    """
    Test init from a m80 object.
    Needs to inherit from m80_Fasta,
    but dont use it directly, use the m80
    API instead
    """
    # grab from mm80
    f = lp.Fasta("ACGT")
    assert len(f["chr1"]) == 100


def test_init_tables(smpl_fasta):
    tables = set(
        [
            x[0]
            for x in smpl_fasta.m80.db.cursor()
            .execute("SELECT name FROM sqlite_master")
            .fetchall()
        ]
    )
    assert "added_order" in tables
    assert "nicknames" in tables
    assert "attributes" in tables


def test_add_chrom(smpl_fasta):
    chrom = lp.Chromosome("added", "U" * 100)
    smpl_fasta.add_chrom(chrom)
    assert len(smpl_fasta["added"]) == 100
    smpl_fasta.del_chrom("added")


def test_add_duplicate_chrom(smpl_fasta):
    chrom = lp.Chromosome("added", "U" * 100)
    if "added" not in smpl_fasta:
        smpl_fasta.add_chrom(chrom)
    with pytest.raises(ValueError):
        smpl_fasta.add_chrom(chrom)
    smpl_fasta.del_chrom("added")


def test_add_duplicate_chrom_with_replace(smpl_fasta):
    if "added" not in smpl_fasta:
        chrom = lp.Chromosome("added", "U" * 100)
        smpl_fasta.add_chrom(chrom)
    chrom = lp.Chromosome("added", "U" * 200)
    smpl_fasta.add_chrom(chrom, replace=True)
    assert len(smpl_fasta["added"]) == 200
    smpl_fasta.del_chrom("added")


def test_remove_chrom(smpl_fasta):
    if "added" not in smpl_fasta:
        chrom = lp.Chromosome("added", "U" * 100)
        smpl_fasta.add_chrom(chrom)
    smpl_fasta.del_chrom("added")
    assert "added" not in smpl_fasta


def test_remove_Chromosome(smpl_fasta):
    if "added" not in smpl_fasta:
        chrom = lp.Chromosome("added", "U" * 100)
        smpl_fasta.add_chrom(chrom)
    chrom = smpl_fasta["added"]
    smpl_fasta.del_chrom(chrom)
    assert "added" not in smpl_fasta


def test_remove_wrong_type(smpl_fasta):
    with pytest.raises(ValueError):
        # must be a str or a Chromosome
        smpl_fasta.del_chrom(4)


def test_remove_chrom_missing(smpl_fasta):
    with pytest.raises(ValueError):
        smpl_fasta.del_chrom("abcdefg")


def test_from_file(m80_Fasta):
    # this fixture creates an m80 Fasta
    # and returns true
    assert m80_Fasta == True


def test_iter(smpl_fasta):
    for c in smpl_fasta:
        assert isinstance(c, lp.Chromosome)


def test__len__(smpl_fasta):
    chroms = list(smpl_fasta)
    assert len(chroms) == len(smpl_fasta)


def test__contains__(smpl_fasta):
    assert "chr1" in smpl_fasta


def test__getitem__(smpl_fasta):
    assert isinstance(smpl_fasta["chr1"], lp.Chromosome)


def test_chrom_attrs(smpl_fasta, m80_Fasta):
    c1 = smpl_fasta["chr1"]
    assert "foobar" in c1._attrs
    c2 = lp.Fasta("ACGT")
    assert "foo" in c2["chr1"]._attrs


def test_add_attribute(smpl_fasta):
    smpl_fasta._add_attribute("chr1", "added")
    c1 = smpl_fasta["chr1"]
    assert "added" in c1._attrs


def test_add_nickname(smpl_fasta):
    smpl_fasta._add_nickname("chr1", "NICKNAME")
    n1 = smpl_fasta["NICKNAME"]
    c1 = smpl_fasta["chr1"]
    assert n1 == c1


def test_get_nickname(smpl_fasta):
    n1 = smpl_fasta["CHR1"]
    c1 = smpl_fasta["chr1"]
    assert n1 == c1


def test_get_chrom_names(smpl_fasta):
    names = smpl_fasta.chrom_names()
    for x in ["chr1", "chr2", "chr3", "chr4"]:
        assert x in names


def test_contains_with_chromosome_object(smpl_fasta):
    chr1 = smpl_fasta["chr1"]
    assert chr1 in smpl_fasta


def test_get_chrom_not_in_fasta(smpl_fasta):
    with pytest.raises(ValueError):
        smpl_fasta["Nope"]


def test_to_Fasta_file(smpl_fasta):
    # this does not test that the fasta file is
    # correct, just that the
    tfile = smpl_fasta.m80.tmpfile()
    smpl_fasta.to_fasta(tfile.name)
    # now read it back into a new Fasta object
    if m80.exists("Fasta", "copy"):
        m80.delete("Fasta", "copy")
    fasta_copy = lp.Fasta.from_file("copy", tfile.name)
    for chrom in smpl_fasta:
        assert chrom in fasta_copy
    m80.delete("Fasta", "copy")
    # Delete the copy
