import pytest
import os
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
            'raw', 'ZmB73_5b_FGS.gff.gz'
        )
    )
    # This is stupid and necessary because pytables wont let me open
    # more than one table
    x = Loci('Zm5bFGS')
    if len(x) == 0:
        x.add_gff(  
            gff
        )
    return x


