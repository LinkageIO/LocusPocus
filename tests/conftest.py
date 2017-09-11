import pytest
import os
from locuspocus import Locus
from locuspocus import Loci

@pytest.fixture(scope='module')
def simpleLoci():
    # Create a Locus
    a = Locus(1,100,150, id='gene_a')
    # Create a couple more!
    b = Locus(1,160,175, id='gene_b')
    c = Locus(1,180,200, id='gene_c')
    d = Locus(1,210,300, id='gene_d')
    e = Locus(2,100,150, id='gene_e')

    x = Loci('simpleLoci')
    x.add_loci([a,b,c,d,e])
    return x

@pytest.fixture(scope="module")
def testRefGen():
    # We have to build it
    gff = os.path.expanduser(
        os.path.join(
            'raw', 'ZmB73_5b_FGS.gff.gz'
        )
    )
    x = Loci('Zm5bFGS')
    if len(x) == 0:
        x.add_gff(  
            gff
        )
    return x


