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
