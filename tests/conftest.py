import pytest
from locuspocus import Locus
from locuspocus import Loci

@pytest.fixture(scope='module')
def testRefGen():
    # Create a Locus
    a = Locus(8,100,150, name='gene_a')
    # Createa  couple more!
    b = Locus(8,160,175, name='gene_b')
    c = Locus(8,180,200, name='gene_c')
    d = Locus(8,210,300, name='gene_d')
    e = Locus(9,100,150, name='gene_e')

    x = Loci()
    x.append([a,b,c,d,e])
    return x

''' -------------------------------------------------------------------------
            'Test' Fixtures
'''

@pytest.fixture(scope='module')
def testRefGen(Zm5bFGS):
    # This was a mistake
    return Zm5bFGS

''' -------------------------------------------------------------------------
            RefGen Fixtures
'''

@pytest.fixture(scope="module")
def testRefGen():
    # We have to build it
    gff = os.path.expanduser(
        os.path.join(
            cf.options.testdir,
            'raw', 'ZmB73_5b_FGS.gff.gz'
        )
    )
    # This is stupid and necessary because pytables wont let me open
    # more than one table
    co.RefGen.from_gff(
        gff, 'Zm5bFGS'
    )
    return co.RefGen('Zm5bFGS')


