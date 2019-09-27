import pytest
import os

from locuspocus import Locus
from locuspocus import FrozenLoci
from locuspocus import Fasta

import minus80.Tools as m80tools

from locuspocus import Chromosome

@pytest.fixture(scope='module')
def simpleFrozenLoci():
    m80tools.delete('FrozenLoci','simpleFrozenLoci')
    # Create a Locus
    a = Locus(1,100,150, id='gene_a')
    # Create a couple more!
    b = Locus(1,160,175, id='gene_b')
    c = Locus(1,180,200, id='gene_c')
    d = Locus(1,210,300, id='gene_d')
    e = Locus(2,100,150, id='gene_e')

    x = FrozenLoci('simpleFrozenLoci')
    x.add_loci([a,b,c,d,e])
    return x

@pytest.fixture(scope="module")
def testFrozenLoci():
    # We have to build it
    if m80tools.available('FrozenLoci','Zm5bFGS'):
        x =  FrozenLoci('Zm5bFGS')
    else:
        m80tools.delete('FrozenLoci','Zm5bFGS')
        gff = os.path.expanduser(
            os.path.join(
                'raw', 
                'ZmB73_5b_FGS.gff.gz'
            )
        )
        x = FrozenLoci('Zm5bFGS')
        if len(x) == 0:
            x.import_gff(  
                gff
            )
    x.set_primary_feature_type('gene')
    return x


@pytest.fixture(scope='module')
def m80_Fasta():
    '''
        Create a Fasta which doesn't get 
        returned. Access the Fasta through
        the m80 API
    '''
    # delete the onl
    m80tools.delete('Fasta','ACGT')
    f = Fasta.from_file('ACGT','raw/ACGT.fasta')
    return True


@pytest.fixture(scope='module')                                                               
def smpl_fasta():                                                                  
    ''' A simple fasta that agrees with smpl_annot'''                           
    m80tools.delete('Fasta','smpl_fasta')
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
