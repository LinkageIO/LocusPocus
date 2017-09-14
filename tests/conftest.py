import pytest
import os
from locuspocus import Locus
from locuspocus import Loci
from locuspocus import Fasta

from locuspocus.Fasta import Chromosome

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

@pytest.fixture()                                                               
def smpl_fasta():                                                                  
    ''' A simple fasta that agrees with smpl_annot'''                           
    fasta = Fasta()                                                             
    chr1 = Chromosome('A'*500000)                                                  
    chr2 = Chromosome('C'*500000)                                               
    chr3 = Chromosome('G'*500000)                                               
    chr4 = Chromosome('T'*500000)                                                  
    fasta.add_chrom('chr1',chr1)                                                
    fasta.add_chrom('chr2',chr2)                                                
    fasta.add_chrom('chr3',chr3)                                                
    fasta.add_chrom('chr4',chr4)                                                
    return fasta  
