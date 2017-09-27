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


def test_to_minus_80(smpl_fasta):
    smpl = smpl_fasta
    smpl.to_minus80('smpleFasta')
    thawed = lp.Fasta.from_minus80('smpleFasta')
    assert smpl.added_order == thawed.added_order
    for chrom,seq in smpl.chroms.items():
        assert chrom in thawed
        assert seq.seq == thawed[chrom].seq
    for k,v in smpl.nicknames.items():
        assert k in thawed.nicknames
        assert thawed.nicknames[k] == v
    for k,v in smpl.attributes.items():
        assert k in thawed.attributes
        assert thawed.attributes[k] == v


       
