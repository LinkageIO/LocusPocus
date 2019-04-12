import pytest
import os

from locuspocus import Locus
from locuspocus import RefLoci
from locuspocus import Fasta

import minus80.Tools as m80tools

from locuspocus.Fasta import Chromosome

@pytest.fixture(scope='module')
def simpleRefLoci():
    m80tools.delete('RefLoci','simpleRefLoci',force=True)
    # Create a Locus
    a = Locus(1,100,150, id='gene_a')
    # Create a couple more!
    b = Locus(1,160,175, id='gene_b')
    c = Locus(1,180,200, id='gene_c')
    d = Locus(1,210,300, id='gene_d')
    e = Locus(2,100,150, id='gene_e')

    x = RefLoci('simpleRefLoci')
    x.add_loci([a,b,c,d,e])
    return x

@pytest.fixture(scope="module")
def testRefGen():
    # We have to build it
    if m80tools.available('RefLoci','Zm5bFGS'):
        return RefLoci('Zm5bFGS')
    m80tools.delete('RefLoci','Zm5bFGS',force=True)
    gff = os.path.expanduser(
        os.path.join(
            'raw', 'ZmB73_5b_FGS.gff.gz'
        )
    )
    x = RefLoci('Zm5bFGS')
    if len(x) == 0:
        x.import_gff(  
            gff
        )
    return x


@pytest.fixture(scope='module')
def m80_Fasta():
    '''
        Create a Fasta which doesn't get 
        returned. Access the Fasta through
        the m80 API
    '''
    # delete the onl
    m80tools.delete('Fasta','ACGT',force=True)
    f = Fasta.from_file('ACGT','raw/ACGT.fasta')
    return True


@pytest.fixture(scope='module')                                                               
def smpl_fasta():                                                                  
    ''' A simple fasta that agrees with smpl_annot'''                           
    m80tools.delete('Fasta','smpl_fasta',force=True)
    fasta = Fasta('smpl_fasta')                                                             
    chr1 = Chromosome('chr1','A'*500000)                                                  
    chr2 = Chromosome('chr2','C'*500000)                                               
    chr3 = Chromosome('chr3','G'*500000)                                               
    chr4 = Chromosome('chr4','T'*500000)                                                  
    fasta.add_chrom(chr1)                                                
    fasta.add_chrom(chr2)                                                
    fasta.add_chrom(chr3)                                                
    fasta.add_chrom(chr4)                                                
    fasta._add_nickname('chr1','CHR1')
    fasta._add_attribute('chr1','foobar')
    return fasta  
