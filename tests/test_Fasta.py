'''
    Tests
'''
import pytest

def test_init(smpl_fasta):
    assert len(smpl_fasta['chr1']) == 500000
    assert len(smpl_fasta['chr2']) == 500000
    assert len(smpl_fasta['chr3']) == 500000
    assert len(smpl_fasta['chr4']) == 500000
