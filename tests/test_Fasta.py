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


#def test_minus80_init(smpl_fasta):
#    smpl = smpl_fasta
#    thawed = lp.Fasta('smpl_fasta')
#    assert list(smpl) == list(thawed)

       
