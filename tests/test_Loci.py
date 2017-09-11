import pytest

from locuspocus import Locus

'''
    Unit tests
'''

def test_get_loci_from_ids(testRefGen):
    random_loci = sorted(testRefGen.random_loci(n=10))
    get_loci_from_ids = sorted(testRefGen.from_ids([x.id for x in random_loci]))
    assert set(random_loci) == set(get_loci_from_ids)

def test_get_item(testRefGen):
    random_locus = testRefGen.random_locus()
    assert random_locus == testRefGen[random_locus.id]

def test_get_items_from_list(testRefGen):
    random_loci = sorted(testRefGen.random_loci(n=10))
    get_loci_from_ids = sorted(testRefGen[[x.id for x in random_loci]])
    assert set(random_loci) == set(get_loci_from_ids)

def test_lowercase_get_item(testRefGen):
    random_locus = testRefGen.random_locus()
    random_id = random_locus.id
    # Stupid mutability
    random_id.lower()
    assert random_locus == testRefGen[random_id]

def test_loci_within(testRefGen):
    random_locus = testRefGen.random_locus()
    bigger_locus = Locus(
        random_locus.chrom,
        start=random_locus.start-100,
        end=random_locus.end+100
    )
    genes = testRefGen.genes_within(bigger_locus)
    assert random_locus in genes

def test_locus_not_in_upstream_downstream(testRefGen):
    '''
        Upstream and downstream should not include the gene of interest.
    '''
    random_locus = testRefGen.random_locus()
    upstream = testRefGen.upstream_loci(
        random_locus,window_size=50e5,locus_limit=5
    )
    downstream = testRefGen.downstream_loci(
        random_locus,locus_limit=5,window_size=50e5
    )
    assert random_locus not in upstream
    assert random_locus not in downstream

def test_upstream_downstream_loci(testRefGen):
    '''
        Take downstream of genes, then upstream genes of the 
        last gene in downstream. Tests that returns the same interval 
    '''
    # Grab downstream genes of random genes
    random_locus = testRefGen.random_locus()
    # Grab 10 downstream genes
    downstream_loci = testRefGen.downstream_loci(
        random_locus,locus_limit=11,window_size=50e10
    )
    assert len(downstream_loci) == 11
    # grab last gene
    last_gene = downstream_loci.pop(-1)
    # Grab upstream genes
    upstream_loci = testRefGen.upstream_loci(
        last_gene,locus_limit=10,window_size=50e10
    )
    assert sorted(downstream_loci) == sorted(upstream_loci)

def test_flanking_loci(testRefGen):
    random_locus = testRefGen.random_locus()
    downstream = testRefGen.downstream_loci(
        random_locus, window_size=50e6, locus_limit=5
    )
    upstream = testRefGen.upstream_loci(
        random_locus, window_size=50e6, locus_limit=5
    )
    flanking = testRefGen.flanking_loci(
        random_locus, window_size=50e6, flank_limit=5
    )
    assert sorted(flanking) == sorted(upstream + downstream)

def test_flanking_loci_includes_within_loci_for_SNPS(testRefGen):
    random_locus = testRefGen.random_locus()
    # test snp
    test_snp = Locus(random_locus.chrom,random_locus.start,window=50e5)
    flanking = testRefGen.flanking_loci(test_snp)
    assert random_locus not in flanking

def test_candidate_loci_from_SNP(testRefGen):
    random_locus = testRefGen.random_locus()
    # grab a bunch of downstream genes
    down1,down2 = testRefGen.downstream_loci(
        random_locus,locus_limit=2,window_size=50e6
    )
    # Create a Locus that is on gene 5
    test_snp = Locus(
        down1.chrom,
        down1.start-50,
        end=down2.end+50,
        window=50e6
    )
    candidates = testRefGen.candidate_loci(
        test_snp,flank_limit=5,chain=False
    )
    assert len(candidates) == 12 

def test_candidate_loci_from_gene_includes_gene(testRefGen):
    random_locus = testRefGen.random_locus()
    # grab a bunch of downstream genes
    downstream = testRefGen.downstream_loci(
        random_locus,locus_limit=10,window_size=50e6
    )
    # Create a Locus that is on gene 5
    candidates = testRefGen.candidate_loci(
        downstream[5],flank_limit=10,window_size=50e6
    )
    assert downstream[4] in candidates

def test_non_chained_candidates(testRefGen):
    random_loci = testRefGen.random_loci(n=10)
    # Create a Locus that is on gene 5
    candidates = testRefGen.candidate_loci(
        random_loci,flank_limit=10,window_size=50e6,chain=False
    )
    # test that we got candidates for each random locus
    assert len(candidates) == len(random_loci)
   

def test_flank_limit_for_candidate_loci(testRefGen):
    random_locus = testRefGen.random_locus()
    # Create a Locus that is on gene 5
    candidates = testRefGen.candidate_loci(
        random_locus,flank_limit=5,window_size=50e6,chain=True
    )
    assert len(candidates) == 11

def test_flank_limit_for_candidate_loci_from_SNP(testRefGen):
    random_locus = testRefGen.random_locus()
    downstream = testRefGen.downstream_loci(
        random_locus,locus_limit=10,window_size=50e6
    )
    test_snp = Locus(downstream[5].chrom,downstream[5].start,window=50e6)
    # Create a Locus that is on gene 5
    candidates = testRefGen.candidate_loci(
        test_snp,flank_limit=5,window_size=50e6
    )
    assert len(candidates) == 11

def test_bootstrap_candidate_length_equal_from_SNP(testRefGen):
    random_locus = testRefGen.random_locus()
    test_snp = Locus(random_locus.chrom,random_locus.start,window=50e6)
    candidates = testRefGen.candidate_loci(test_snp)
    bootstraps = testRefGen.bootstrap_candidate_loci(test_snp)
    assert len(candidates) == len(bootstraps)

def test_bootstrap_candidate_length_equal_from_gene(testRefGen):
    random_locus = testRefGen.random_locus()
    candidates = testRefGen.candidate_loci(random_locus,window_size=5e10)
    bootstraps = testRefGen.bootstrap_candidate_loci(random_locus,window_size=5e10)
    assert len(candidates) == len(bootstraps)

def test_refgen_length(testRefGen):
    # grab length from sqlite 
    from_sql = testRefGen._db.cursor().execute('''
        SELECT COUNT(*) FROM loci;
    ''').fetchone()[0]
    assert from_sql == len(testRefGen)

def test_random_loci_returns_correct_n(testRefGen):
    assert len(testRefGen.random_loci(n=50)) == 50

# New Tests

def test_add_loci(simpleLoci):
    new_locus = Locus(1,100,300,id='new_locus')
    simpleLoci.add_loci(new_locus)
    assert new_locus in simpleLoci

def test_removie_locus(simpleLoci):
    new_locus = Locus(1,100,300,id='new_locus')
    if new_locus not in simpleLoci:
        simpleLoci.add_locus(new_locus)
    assert new_locus in simpleLoci
    simpleLoci.remove_locus(new_locus.id)
    assert new_locus not in simpleLoci

def test_add_gff(testRefGen):
    assert len(testRefGen) > 0

def test_len(simpleLoci):
    assert isinstance(len(simpleLoci),int)

def test_num_loci(simpleLoci):
    assert isinstance(simpleLoci.num_loci(),int)

def test_random_locus(simpleLoci):
    assert isinstance(simpleLoci.random_locus(),Locus)

def test_random_loci(simpleLoci):
    loci = simpleLoci.random_loci(n=3)
    assert len(loci) == 3
    for x in loci:
        assert isinstance(x,Locus)

def test_iter_loci(simpleLoci):
    for x in simpleLoci.iter_loci():
        assert isinstance(x,Locus)

def test_intersection():
    pass

def test_from_id(simpleLoci):
    locus = simpleLoci.get_locus_from_id('gene_a')
    assert locus.id in simpleLoci
    assert isinstance(locus,Locus)

def test_get_loci_from_ids(simpleLoci):
    loci = simpleLoci.get_loci_from_ids(['gene_a','gene_b','gene_c'])
    assert len(loci) == 3

def test_get_item(simpleLoci):
    assert isinstance(simpleLoci['gene_a'],Locus)

def test_encompassing_loci(simpleLoci):
    pass

def test_loci_within(simpleLoci):
    pass

def test_upstream_loci(simpleLoci):
    pass

def test_downstream_loci(simpleLoci):
    pass

def test_flanking_loci(simpleLoci):
    pass

def test_candidate_loci(simpleLoci):
    pass

def test_bootstrap_candidate_loci(simpleLoci):
    pass

def test_pairwise_distance(simpleLoci):
    pass

def test_contains(simpleLoci):
    assert 'gene_a' in simpleLoci 

def test_add_alias(simpleLoci):
    simpleLoci._remove_aliases()
    assert 'alias_1' not in simpleLoci
    simpleLoci.add_alias('gene_a','alias_1')
    assert 'alias_1' in simpleLoci

def test_num_aliases(simpleLoci):
    simpleLoci._remove_aliases()
    simpleLoci.add_alias('gene_a','alias_1')
    assert simpleLoci.num_aliases() == 1
    simpleLoci._remove_aliases()

def test_aliases(simpleLoci):
    simpleLoci._remove_aliases()
    simpleLoci.add_alias('gene_a','alias_1')
    x = simpleLoci['alias_1']
    assert x.id == 'gene_a'
    simpleLoci._remove_aliases()

def test_remove_aliases(simpleLoci):
    simpleLoci._remove_aliases()
    assert simpleLoci.num_aliases() == 0

def test_has_annotations():
    pass

def test_export_annotations():
    pass

def test_add_annotations():
    pass

def test_remove_annotations():
    pass

