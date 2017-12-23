'''
    Tests
'''
import pytest
import locuspocus as lp

def test_init(smpl_fasta):
    assert len(smpl_fasta['chr1']) == 500000
    assert len(smpl_fasta['chr2']) == 500000
    assert len(smpl_fasta['chr3']) == 500000
    assert len(smpl_fasta['chr4']) == 500000


def test_m80_init(m80_Fasta):
    '''
    Test init from a m80 object.
    Needs to inherit from m80_Fasta,
    but dont use it directly, use the m80
    API instead
    '''
    # grab from mm80
    f = lp.Fasta('ACGT')
    assert len(f['chr1']) == 100
      
def test_init_tables(smpl_fasta):
    tables = set([x[0] for x in smpl_fasta._db.cursor().execute(
        'SELECT name FROM sqlite_master'        
    ).fetchall()])
    assert 'added_order' in tables
    assert 'nicknames' in tables
    assert 'attributes' in tables

def test_add_chrom(smpl_fasta):
    chrom = lp.Chromosome('added','U'*100)
    smpl_fasta.add_chrom(chrom)
    assert len(smpl_fasta['added']) == 100

def test_from_file(m80_Fasta):
    # this fixture creates an m80 Fasta
    # and returns true
    assert m80_Fasta == True

def test_iter(smpl_fasta):
    for c in smpl_fasta:
        assert isinstance(c,lp.Chromosome)

def test__len__(smpl_fasta):
    chroms = list(smpl_fasta)
    assert len(chroms) == len(smpl_fasta)

def test__contains__(smpl_fasta):
    assert 'chr1' in smpl_fasta 

def test__getitem__(smpl_fasta):
    assert isinstance(smpl_fasta['chr1'],lp.Chromosome)

def test_chrom_attrs(smpl_fasta,m80_Fasta):
    c1 = smpl_fasta['chr1']
    assert 'foobar' in c1._attrs
    c2 = lp.Fasta('ACGT')
    assert 'foo' in c2['chr1']._attrs

def test_add_attribute(smpl_fasta):
    smpl_fasta._add_attribute('chr1','added')
    c1 = smpl_fasta['chr1']
    assert 'added' in c1._attrs

def test_add_nickname(smpl_fasta):
    smpl_fasta._add_nickname('chr1','NICKNAME')
    n1 = smpl_fasta['NICKNAME']
    c1 = smpl_fasta['chr1']
    assert n1 == c1

def test_get_nickname(smpl_fasta):
    n1 = smpl_fasta['CHR1']
    c1 = smpl_fasta['chr1']
    assert n1 == c1
    
