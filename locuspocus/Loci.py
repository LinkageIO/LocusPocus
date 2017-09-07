#!/usr/bin/python3
import itertools
import collections
import random
import numpy as np
import scipy as sp
import math
import gzip
import re
import logging 

from minus80 import Freezable
from collections import defaultdict
from functools import lru_cache,wraps

from .LocusDist import LocusDist
from .Locus import Locus
from .Exceptions import ZeroWindowError

def memoize(obj):
    cache = obj.cache = {}
    @wraps(obj)
    def memoizer(*args, **kwargs):
        # Give us a way to clear the cache
        if 'clear_cache' in kwargs:
            cache.clear()
        # This wraps the calling of the memoized object
        key = str(args) + str(kwargs)
        if key not in cache:
            cache[key] = obj(*args, **kwargs)
        return cache[key]
    return memoizer


class Loci(Freezable):
    '''
        Just a bunch of Locuses. Loci are more than the sum of their 
        parts. They have a name and represent something bigger than
        theirselves. They are important. They live on the disk in a
        database.
    '''
    log = logging.getLogger(__name__)
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
                    '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(logging.INFO)

    def __init__(self,name):
        super().__init__(name)
        self._initialize_tables()

    def add_loci(self,loci):
        if isinstance(loci,Locus):
            # Its actually just one locus
            locus = loci
            self._db.cursor().execute('''
            INSERT OR REPLACE INTO loci VALUES (?,?,?,?)
            ''',(locus.name, locus.chrom, locus.start, locus.end))
            self._db.cursor().executemany('''
            INSERT OR REPLACE INTO loci_attrs VALUES (?,?,?)
            ''',[(locus.id,key,val) for key,val in locus.attr.items()])
        else:
            try:
                # support adding lists of genes
                self.log.info('Adding {} Genes info to database'.format(len(loci)))
                cur = self._db.cursor()
                cur.execute('BEGIN TRANSACTION')
                cur.executemany(
                    'INSERT OR REPLACE INTO loci VALUES (?,?,?,?)',
                    ((x.name,x.chrom,x.start,x.end) for x in loci)
                )
                self.log.info('Adding Gene attr info to database')
                cur.executemany(
                    'INSERT OR REPLACE INTO loci_attrs VALUES (?,?,?)',
                    ((x.id,key,val) for x in loci for key,val in x.attr.items())
                )
                cur.execute('END TRANSACTION')
            except Exception as e:
                cur.exectue('ROLLBACK')

    def add_gff(self,filename,
                 chrom_feature='chromosome',gene_feature='gene',
                 ID_attr='ID',attr_split='='):
        '''
            Imports Loci from a gff (General Feature Format) file.
            See more about the format here:
            http://www.ensembl.org/info/website/upload/gff.html

            Parameters
            ----------

            filename : str
                The path to the GFF file.
            name : str
                The name if the RefGen object to be stored in the core
                camoco database.
            description : str
                A short description of the RefGen for future reference
            chrom_feature : str (default: chromosome)
                The name of the feature (in column 3) that designates a
                a chromosome.
            gene_feature : str (default: gene)
                The name of the feature (in column 2) that designates a 
                gene. These features will be the main object that the RefGen
                encompasses. 
           [ID_attr : str (default: ID)
                The key in the attribute column which designates the ID or 
                name of the feature.
            attr_split : str (default: '=')
                The delimiter for keys and values in the attribute column
        '''
        loci = list()
        chroms = dict()
        if filename.endswith('.gz'):
            IN = gzip.open(filename,'rt')
        else:
            IN = open(filename,'r')
        for line in IN:
            #skip comment lines
            if line.startswith('#'):
                continue
            (chrom,source,feature,start,
             end,score,strand,frame,attributes) = line.strip().split('\t')
            attributes = dict([(field.strip().split(attr_split)) \
                    for field in attributes.strip(';').split(';')])
            if feature == gene_feature:
                loci.append(
                    Locus(
                        chrom,int(start),int(end),
                        attributes[ID_attr].upper().strip('"'),strand=strand,
                        **attributes
                    ).update(attributes)
                )
                # Keep track of seen chromosomes
                if chrom not in chroms:
                    chroms[chrom] = end
                else:
                    if end > chroms[chrom]:
                        chroms[chrom] = end
        IN.close()
        self.add_loci(loci)


    ###########################################
    ########### Internal Methods  #############
    ###########################################

    def __len__(self):
        '''
            Returns the number of Loci in the database
        '''
        return self.num_loci()

    def _initialize_tables(self):
        '''
            Initializes the Tables holding all the information
            about the Loci. 
        '''
        cur = self._db.cursor()
        cur.execute('''
            CREATE TABLE IF NOT EXISTS loci (
                id TEXT NOT NULL UNIQUE,
                chromosome TEXT NOT NULL,
                start INTEGER,
                end INTEGER
            );
            /*
            ----Create a table that contains gene attribute
                mappings
            */
            CREATE TABLE IF NOT EXISTS loci_attrs (
                id TEXT NOT NULL,
                key TEXT,
                val TEXT
            );
            CREATE TABLE IF NOT EXISTS aliases (
                alias TEXT UNIQUE,
                id TEXT
            );
            CREATE TABLE IF NOT EXISTS func (
                id TEXT,
                desc TEXT,
                UNIQUE(id,desc) ON CONFLICT IGNORE
            );
            CREATE TABLE IF NOT EXISTS ortho_func (
                id TEXT,
                desc TEXT,
                UNIQUE(id,desc) ON CONFLICT IGNORE
            );
            
            CREATE INDEX IF NOT EXISTS gene_start_end ON loci (chromosome,start DESC, end ASC, id);
            CREATE INDEX IF NOT EXISTS gene_end_start ON loci (chromosome,end DESC,start DESC,id);
            CREATE INDEX IF NOT EXISTS gene_start ON loci (chromosome,start);
            CREATE INDEX IF NOT EXISTS gene_end ON loci (chromosome,end);
            CREATE INDEX IF NOT EXISTS geneid ON loci (id);
            CREATE INDEX IF NOT EXISTS geneattr ON loci_attrs (id);
            CREATE INDEX IF NOT EXISTS id ON func(id);
            CREATE INDEX IF NOT EXISTS id ON ortho_func(id); 
            
            ''');

    @memoize
    def num_loci(self):
        '''
            Returns the number of genes in the dataset
        '''
        return self._db.cursor().execute(
           ''' SELECT COUNT(*) FROM loci'''
        ).fetchone()[0]

    def random_locus(self,**kwargs):
        '''
            Returns a random locus within loci.
            Also allows passing of keyword arguments to Locus
            constructor method allowing for flexible generation.
            See Locus.__init__ for more details.

            Parameters
            ----------
            **kwargs : key,value pairs
                Extra parameters passed onto the locus init method.

            Returns
            -------
            A Locus object (camoco.Locus based)

        '''
        random_rowid = random.randint(1,self.num_loci())
        chrom,start,end,id = self._db.cursor().execute('''
            SELECT chromosome,start,end,id from loci WHERE rowid = ?
        ''',(random_rowid,)).fetchone()
        return Locus(chrom,start,end,id=id,**kwargs
        )


    def intersection(self,loci):
        '''
            Return the subset of loci that are in the dataset.
            
            Parameters
            ----------
            loci : list-like of Locus

            Returns
            -------
            a list like object containing loci that
            are in the dataset.
        '''
        return [x for x in loci if x in self]

    def random_loci(self,n,**kwargs):
        '''
            Return random loci from the database, without replacement. 

            Parameters
            ----------
            n : int

            **kwargs : key,value pairs
                Extra parameters passed onto the locus init method

            Returns
            -------
            An iterable containing n (unique) random loci

        '''
        rand_nums = np.random.choice(self.num_loci()+1,n,replace=False)
        gene_info = self._db.cursor().executemany(
                "SELECT chromosome,start,end,id from loci WHERE rowid = ?",
                [[int(rownum)] for rownum in rand_nums]
        )
        return set([Locus(chr,start,end=end,id=id,**kwargs) for \
            (chr,start,end,id) in gene_info])

    ## UnChecked

    def iter_genes(self):
        '''
            Iterates over genes in RefGen.

            Returns
            -------
            A generator containing genes
        '''
        for x in self._db.cursor().execute('''
                SELECT chromosome,start,end,id FROM loci
            '''):
            yield Locus(*x)


    @lru_cache(maxsize=131072)
    def from_id(self, gene_id):
        '''
            Returns a gene object from a string

            Parameters
            ----------
            gene_id : str
                ID of the gene you want to pull out

            Returns
            -------
            A single gene Object

        '''
        cur = self._db.cursor()
        if gene_id not in self:
            result = cur.execute('SELECT id FROM aliases WHERE alias = ?', [gene_id]).fetchone()
            if not result:
                raise ValueError('{} not in {}'.format(gene_id,self.name))
            gene_id = result[0]
        gene_id = gene_id.upper()
        info = cur.execute('SELECT chromosome,start,end,id FROM loci WHERE id = ?', [gene_id]).fetchone()
        return Locus(*info)

    def from_ids(self, gene_ids, check_shape=False):
        '''
            Returns a list of gene object from an iterable of id strings
            OR from a single gene id string.

            Parameters
            ----------
            gene_list : str OR iterable of str
                ID(s) of the genes you want to pull out
            check_shape : bool (default: False)
                Check if you get back the same number of ids you
                pass in. If false (default), just give back what
                you find, ignoring erronous ids.

            Returns
            -------
            A list of locus objects if you pass in an iterable,
            otherwise a single gene

        '''
        if isinstance(gene_ids,str):
            import warnings
            warnings.warn(
                'Passing singe values into RefGen.from_ids is deprecated. Use RefGen.from_id() '
                'or slicing syntax instead.'    
            )
            return self.from_id(gene_ids)
        genes = []
        for id in gene_ids:
            try:
                genes.append(self.from_id(id))
            except ValueError as e:
                if check_shape == False:
                    continue
                else:
                    raise e
        return genes

    # NOTE: Dont LRU cache this, it gets cached in from_id
    def __getitem__(self,item):
        '''
            A convenience method to extract loci from the reference genome.
        '''
        if isinstance(item,str):
            return self.from_id(item)
        # Allow for iterables of items to be passed in
        else:
            return self.from_ids(item)

    def encompassing_genes(self,loci,chain=True):
        '''
            Returns the gene encompassing the locus. In other words
            if a locus (e.g. a SNP) is inside of a gene, i.e. the
            start of the gene is upstream and the end of the gene
            is downstream of the locus boundaries, this method will
            return it. Not to be confused with candidate genes, which
            will return genes upstream and downstream surrounding a locus.

            Parameters
            ----------
            loci : an iterable of loci
                The method will return encompassing genes for each
                locus in ther iterable. If a single locus is passed,
                a single result will be returned.
            chain : bool (defualt: True)
                Specifies whether or not to chain results. If true
                a single list will be returned, if False, a result for
                each locus passed in will be returned.
        '''
        if isinstance(loci,Locus):
            return [
                Locus(*x) \
                for x in self._db.cursor().execute('''
                    SELECT chromosome,start,end,id FROM loci
                    WHERE chromosome = ?
                    AND start <= ? AND end >= ?
                    ''',
                    (loci.chrom,loci.start,loci.end))
            ]
        else:
            iterator = iter(loci)
            genes = [self.encompassing_genes(locus,chain=chain) for locus in iterator]
            if chain:
                genes = list(itertools.chain(*genes))
            return genes

    def genes_within(self,loci,chain=True):
        '''
            Returns the genes that START within a locus 
            start/end boundry.

            Looks like: (y=yes,returned; n=no,not returned)

            nnn  nnnnnnn  yyyyyy   yyyyy  yyyyyy yyyyyyy
                    start                        end
                -----x****************************x-------------

        '''
        if isinstance(loci,Locus):
            return [
                Locus(*x) \
                for x in self._db.cursor().execute('''
                    SELECT chromosome,start,end,id FROM loci
                    WHERE chromosome = ?
                    AND start >= ? AND start <= ?
                    ''',
                    (loci.chrom,loci.start,loci.end))
            ]
        else:
            iterator = iter(loci)
            genes = [self.genes_within(locus,chain=chain) for locus in iterator]
            if chain:
                genes = list(itertools.chain(*genes))
            return genes

    def upstream_genes(self,locus,gene_limit=1000,window_size=None):
        '''
            Find genes that START upstream of a locus. 
            Loci are ordered so that the nearest genes are 
            at the beginning of the list.

            Return Genes that overlap with the upstream window,
            This includes partially overlapping genes, but NOT
            genes that are returned by the genes_within method. 

            Looks like: (y=yes,returned; n=no,not returned)

            nnn  yyyyyyy   yyyyyy   yyyyy  yyyyyy nnnn nnnn nnnnnnnn
                                             nnnn
                                           start             end
                -----------------------------x****************x--
                   ^_________________________| Window (upstream)
        '''
        if locus.window == 0 and window_size is None:
            raise CamocoZeroWindowError(
                'Asking for upstream genes for {}',
                locus.id
            )
        if window_size is not None:
            upstream = locus.start - window_size
        else:
            upstream = locus.upstream
        return [
            Locus(*x) \
            for x in self._db.cursor().execute('''
                SELECT chromosome,start,end,id FROM loci 
                INDEXED BY gene_start_end
                WHERE chromosome = ?
                AND start >= ?  -- Gene must end AFTER locus window (upstream) 
                AND start < ? -- Gene must start BEFORE locus
                ORDER BY start DESC
                LIMIT ?
            ''',(locus.chrom, upstream,locus.start, gene_limit)
        )]

    def downstream_genes(self,locus,gene_limit=1000,window_size=None):
        '''
            Returns genes downstream of a locus. Genes are ordered 
            so that the nearest genes are at the beginning of the list.

            Return Genes that overlap with the downstream window,
            This includes partially overlapping genes, but NOT
            genes that are returned by the genes_within method. 

            Looks like: (y=yes,returned; n=no,not returned)

            nnn  nnnnnnn   nnnnnn nnnn  yyyy  yyyyyy yyyy yyyyyy  nnnnn
               start             end
              ---x****************x--------------------------------
                                  |_______________________^ Window (downstream)
        '''
        if locus.window == 0 and window_size is None:
            raise CamocoZeroWindowError(
                'Asking for upstream genes for {} with no window size',
                locus.id
            )
        if window_size is not None:
            downstream = locus.end + window_size
        else:
            downstream = locus.downstream

        return [
            Locus(*x) \
            for x in self._db.cursor().execute('''
                SELECT chromosome,start,end,id FROM loci
                INDEXED BY gene_start_end
                WHERE chromosome = ?
                AND start > ?
                AND start <= ?
                ORDER BY start ASC
                LIMIT ?
            ''',(locus.chrom, locus.end, downstream, gene_limit)
        )]

    def flanking_genes(self, loci, flank_limit=2,chain=True,window_size=None):
        '''
            Returns genes upstream and downstream from a locus
            ** done NOT include genes within locus **
        '''
        if isinstance(loci,Locus):
            # If we cant iterate, we have a single locus
            locus = loci
            if locus.window == 0 and window_size is None:
                raise CamocoZeroWindowError(
                    'Asking for upstream genes for {} and no window size.',
                    locus.id
                )
            upstream_gene_limit = int(flank_limit)
            downstream_gene_limit = int(flank_limit)
            up_genes = self.upstream_genes(
                locus, gene_limit=upstream_gene_limit, window_size=window_size
            )
            down_genes = self.downstream_genes(
                locus, gene_limit=downstream_gene_limit, window_size=window_size
            )
            if chain:
                return list(itertools.chain(up_genes,down_genes))
            return (up_genes,down_genes)
        else:
            iterator = iter(loci)
            genes = [
                self.flanking_genes(locus,flank_limit=flank_limit,window_size=window_size)\
                for locus in iterator
            ]
            if chain:
                genes = list(itertools.chain(*genes))
            return genes

    def candidate_genes(self, loci, flank_limit=2,
        chain=True, window_size=None, include_parent_locus=False,
        include_parent_attrs=False, include_num_intervening=False, 
        include_rank_intervening=False, include_num_siblings=False,
        include_SNP_distance=False,attrs=None,return_table=False):
        '''
            Locus to Gene mapping.
            Return Genes between locus start and stop, plus additional
            flanking genes (up to flank_limit)

            Parameters
            ----------
            loci : camoco.Locus (also handles an iterable containing Loci)
                a camoco locus or iterable of loci
            flank_limit : int (default : 2)
                The total number of flanking genes **on each side**
                considered a candidate surrounding a locus
            chain : bool (default : true)
                Calls itertools chain on results before returning
            window_size : int (default: None)
                Optional parameter used to extend or shorten a locus
                window from which to choose candidates from. If None,
                the function will resort to what is available in the
                window attribute of the Locus.
            include_parent_locus : bool (default: False)
                Optional parameter which will update candidate genes
                'attr' attribute with the id of the parent locus
                which contains it.
            include_parent_attrs : iterable (default: False)
                Optional parameter to include attributes from the parent
                locus. Parent locus attrs specified here will be included. 
                If effective loci is > 1, the maximum value will be
                included. E.g. - including the SNP effect size with 
                candidate genes.
            include_num_intervening : bool (default: False)
                Optional argument which adds an attribute to each 
                candidate genes containing the rank of each gene
                as a function of distance away from the parent 
                locus. (i.e. the closest candidate is 1 and the 
                furthest candidate is n)
            include_rank_intervening : bool (default: False)
                Optional argument which adds the rank of each
                candidatea as a funtion of distance from the parent
                Locus. i.e. The closest gene is ranked 1 and the furthest
                gene is ranked n.
            include_num_siblings : bool (default: False)
                Optional argument which adds an attribute to each
                candidate gene containing the number of total 
                candidates (siblings) identifies at the locus.
            include_SNP_distance : bool (default:False)
                Include the distance from the canadidate gene and
                the parent SNP
            attrs : dict (default: None)
                An optional dictionary which will be updated to each
                candidate genes attr value.
            return_table : bool(default: False)
                If True, return a Pandas table (DataFrame)

            Returns
            -------
            a list of candidate genes (or list of lists if chain is False)

        '''
        if isinstance(loci,Locus):
            # If not an iterator, its a single locus
            locus = loci
            genes_within = self.genes_within(locus)
            up_genes,down_genes = self.flanking_genes(
                locus, flank_limit=flank_limit, chain=False,
                window_size=window_size
            )

            # This always returns candidates together, if 
            # you want specific up,within and down genes
            # use the specific methods
            genes = sorted(itertools.chain(up_genes,genes_within,down_genes))
            #include the number of effective loci
            if include_rank_intervening == True:
                ranks = sp.stats.rankdata([abs(x.center_distance(locus)) for x in genes])
            # Iterate through candidate genes and propagate the 
            # parental info
            for i,gene in enumerate(genes):
                #gene.update({'num_effective_loci':len(locus.sub_loci)})
                # include parent locus id if thats specified
                if include_parent_locus == True:
                    gene.update({'parent_locus':locus.id})
                if include_rank_intervening == True:
                    gene.update({'intervening_rank':ranks[i]})
                # update all the parent_attrs
                if include_parent_attrs and len(include_parent_attrs) > 0:
                    if 'all' in include_parent_attrs: 
                        include_parent_attrs = locus.attr.keys()
                    for attr in include_parent_attrs:
                        attr_name = 'parent_{}'.format(attr)
                        gene.update({attr_name: locus[attr]})
            if include_num_intervening == True:
                num_down = 0
                num_up = 0
                # Sort the genes by their distance from the locus
                genes_with_distances = [(gene,abs(gene.center_distance(locus))) for gene in genes]
                genes_with_distances = sorted(genes_with_distances,key=lambda x: x[1])
                for gene,distance in genes_with_distances:
                    if locus.within(gene):
                        gene.update({'num_intervening':-1})
                    elif gene.center >= locus.center:
                        gene.update({'num_intervening':num_down})
                        num_down += 1
                    elif gene.center <= locus.center:
                        gene.update({'num_intervening':num_up})
                        num_up += 1
            if include_num_siblings == True:
                for gene in genes:
                    gene.update({'num_siblings':len(genes)})
            if include_SNP_distance == True:
                for gene in genes:
                    distance = abs(gene.center_distance(locus))
                    gene.update({'SNP_distance':distance})
            if attrs is not None:
                for gene in genes:
                    gene.update(attrs)
            if return_table == True:
                genes = pd.DataFrame([x.as_dict() for x in genes])
            return genes

        else:
            iterator = iter(sorted(loci))
            genes = [
                # This is becoming a pain in the ass
                self.candidate_genes(
                    locus, 
                    flank_limit=flank_limit,
                    chain=chain,
                    window_size=window_size,
                    include_parent_locus=include_parent_locus,
                    include_parent_attrs=include_parent_attrs,
                    include_num_intervening=include_num_intervening,
                    include_rank_intervening=include_rank_intervening,
                    include_num_siblings=include_num_siblings,
                    include_SNP_distance=include_SNP_distance,
                    return_table=return_table,
                    attrs=attrs
                ) for locus in iterator
            ]
            if chain:
                if return_table:
                    genes = pd.concat(genes)
                else:
                    genes = list(set(itertools.chain(*genes)))
            return genes

    def bootstrap_candidate_genes(self, loci, flank_limit=2,
        chain=True, window_size=None, include_parent_locus=False):
        '''
            Returns candidate genes which are random, but conserves
            total number of overall genes.

            Parameters
            ----------
            loci : camoco.Locus (also handles an iterable containing Loci)
                a camoco locus or iterable of loci
            flank_limit : int (default : 2)
                The total number of flanking genes **on each side**
                considered a candidate surrounding a locus
            chain : bool (default : true)
                Calls itertools chain on results before returning,
            include_parent_locus : bool (default: False)
                Optional parameter which will update candidate genes
                'attr' attribute with the id of the parent locus
                which contains it.

            Returns
            -------
            a list of candidate genes (or list of lists if chain is False)

        '''
        if isinstance(loci,Locus):
            # We now have a single locus
            locus = loci
            # grab the actual candidate genes
            num_candidates = len(
                self.candidate_genes(
                    locus, flank_limit=flank_limit,
                    chain=True, window_size=window_size
                )
            )
            if num_candidates == 0:
                return []
            # Snps a random genes from the genome
            random_gene = self.random_locus()
            # Snag the same number of candidates
            random_candidates = self.upstream_genes(
                random_gene, 
                gene_limit=num_candidates,
                window_size=10e100
            )
            if len(random_candidates) != num_candidates:
                # somehow we hit the end of a chromosome
                # or something, just recurse
                random_candidates = self.bootstrap_candidate_genes(
                    locus,flank_limit=flank_limit,chain=True
                )
            if include_parent_locus == True:
                for gene in random_candidates:
                    gene.update({'parent_locus':random_gene.id})
            return random_candidates
        else:
            # Sort the loci so we can collapse down
            locus_list = sorted(loci)
            seen = set()
            bootstraps = list()
            target = self.candidate_genes(
                locus_list,flank_limit=flank_limit,
                chain=False,window_size=window_size
            )
            target_accumulator = 0
            candidate_accumulator = 0
            #self.log.info('target: {}, loci: {}',len(target),len(locus_list))
            for i,(locus,targ) in enumerate(zip(locus_list,target)):
                # compare downstream of last locus to current locus
                candidates = self.bootstrap_candidate_genes(
                    locus, flank_limit=flank_limit, 
                    chain=True, window_size=window_size,
                    include_parent_locus=include_parent_locus
                )
                # If genes randomly overlap, resample
                while len(seen.intersection(candidates)) > 0:
                    candidates = self.bootstrap_candidate_genes(
                        locus, flank_limit=flank_limit,
                        window_size=window_size, chain=True,
                        include_parent_locus=include_parent_locus
                    )
                # Add all new bootstrapped genes to the seen list
                seen |= set(candidates)
                bootstraps.append(candidates)
            if chain:
                bootstraps = list(seen)
            #self.log.info("Found {} bootstraps",len(bootstraps))
            return bootstraps


    def pairwise_distance(self, gene_list=None):
        '''
            returns a vector containing the pairwise distances between genes
            in gene_list in vector form. See np.squareform for matrix
            conversion.
        '''
        if gene_list is None:
            gene_list = list(self.iter_genes())
        query = '''
                SELECT genes.id, chrom.rowid, start, end FROM loci
                LEFT JOIN chromosomes chrom ON genes.chromosome = chrom.id
                WHERE genes.id in ("{}")
                ORDER BY genes.id
        '''.format('","'.join([g.id for g in gene_list]))
        # extract chromosome row ids and gene start positions for each gene
        positions = pd.DataFrame(
            # Grab the chromosomes rowid because its numeric
            self._db.cursor().execute(query).fetchall(),
            columns=['gene','chrom','start','end']
        ).sort_values(by='gene')
        # chromosome needs to be floats
        positions.chrom = positions.chrom.astype('float')
        # Do a couple of checks
        assert len(positions) == len(gene_list), \
            'Some genes in dataset not if RefGen'
        assert all(positions.gene == [g.id for g in gene_list]), \
            'Genes are not in the correct order!'
        distances = RefGenDist.gene_distances(
            positions.chrom.values,
            positions.start.values,
            positions.end.values
        )
        return distances

    def summary(self):
        print ("\n".join([
            'Reference Genome: {} ',
            '{} genes',
            'Genome:',
            '{}']).format(
                self.name,self.num_loci(),
                self.genome
            )
        )

    def plot_loci(self,loci,filename,flank_limit=2):
        '''
            Plots the loci, windows and candidate genes

            Parameters
            ----------
            loci : iterable of co.Loci
                The loci to print
            filename : str
                The output filename
        '''
        plt.clf()
        # Each chromosome gets a plot
        chroms = set([x.chrom for x in loci])
        # Create a figure with a subplot for each chromosome 
        f, axes = plt.subplots(len(chroms),figsize=(10,4*len(chroms)))
        if len(chroms) == 1:
            axes = [axes]
        # Split loci by chromosome
        chromloci = defaultdict(list)
        for locus in sorted(loci):
            chromloci[locus.chrom].append(locus)
        # iterate over Loci
        seen_chroms = set([loci[0].chrom])
        voffset = 1 # Vertical Offset
        hoffset = 0 # Horizonatal Offset
        current_chrom = 0
        for i,locus in enumerate(loci):
            # Reset the temp variables in necessary
            if locus.chrom not in seen_chroms:
                seen_chroms.add(locus.chrom)
                current_chrom += 1
                voffset = 1
                hoffset = 0
            # Do the access things
            cax = axes[current_chrom]
            cax.set_ylabel('Chrom: '+ locus.chrom)
            cax.set_xlabel('Loci')
            cax.get_yaxis().set_ticks([])
            #cax.get_xaxis().set_ticks([])
            # shortcut for current axis
            #cax.hold(True)
            # place marker for start window
            cax.scatter(hoffset,voffset,marker='>')
            # place marker for start snp
            cax.scatter(hoffset+locus.window,voffset,marker='.',color='blue')
            # Place a marker for middle 
            cax.scatter(hoffset+locus.window+len(locus)/2,voffset,marker='|',color='blue')
            # place marker for stop snp
            cax.scatter(hoffset+locus.window+len(locus),voffset,marker='.',color='blue')
            # place marker for stop snp
            cax.scatter(hoffset+locus.window+len(locus)+locus.window,voffset,marker='<')

            #place markers for sub snps
            for subsnp in locus.sub_loci:
                cax.scatter(
                    hoffset + subsnp.start - locus.start + locus.window,
                    voffset,
                    marker='.',
                    color='blue'
                )

            # place a block for interlocal distance
            cax.barh(
                bottom=voffset,
                width=50,
                height=1,
                left=hoffset+locus.window+len(locus)+locus.window,
                color='red'
            )
            # grab the candidate genes
            for gene in self.candidate_genes(locus,flank_limit=flank_limit):
                cax.barh(
                    bottom=voffset,
                    width = len(gene),
                    height= 1,
                    left=gene.start-locus.start+locus.window,
                    color='red'
                )
            voffset += 5

        plt.savefig(filename)
        return f


    def __contains__(self,obj):
        ''' flexible on what you pass into the 'in' function '''
        if isinstance(obj,Locus):
            # you can pass in a gene object (this expression
            # should ALWAYS be true if you
            # created gene object from this RefGen)
            if self._db.cursor().execute(
                '''SELECT COUNT(*) FROM loci WHERE id = ?''',
                (obj.id.upper(),)).fetchone()[0] == 1:
                return True
            else:
                return False
        elif isinstance(obj,str):
            # Can be a string object
            if self._db.cursor().execute('''
                SELECT COUNT(*) FROM loci WHERE id = ?''',
                (str(obj).upper(),)).fetchone()[0] == 1:
                return True
            else:
                return False
        else:
            raise TypeError('Cannot test for containment for {}'.format(obj))

    def add_aliases(self, alias_file, id_col=0, alias_col=1, headers=True):
        ''' 
            Add alias map to the RefGen 
            
            Parameters
            ----------
            alias_file : string (path)
                The path to the alias file
            id_col : int (default: 0)
                The column containing the gene identifier
            alias_col : int (default: 1)
                The columns containing the alias
            header : bool (default: True)
                A switch stating if there is a header row to ignore or not
        '''
        with rawFile(alias_file) as IN:
            if headers:
                garb = IN.readline()
    
            aliases = []
            self.log.info('Importing aliases from: {}',alias_file)
            for line in IN:
                row = re.split(',|\t',line)
                if row[id_col].strip() in self:
                    aliases.append(
                        (row[alias_col],row[id_col].strip())
                    )
        cur = self._db.cursor()
        self.log.info('Saving them in the alias table.')
        cur.execute('BEGIN TRANSACTION')
        cur.executemany(
            'INSERT OR REPLACE INTO aliases VALUES (?,?)',
            aliases
        )
        cur.execute('END TRANSACTION')
    def num_aliases(self):
        '''
            Returns the number of aliases currently in the database
        '''
        return self._db.cursor() \
            .execute('SELECT COUNT(*) FROM aliases') \
            .fetchone()[0]

    def aliases(self, gene_id):
        if isinstance(gene_id,str):
            return [alias[0] for alias in self._db.cursor().execute('''
                    SELECT alias FROM aliases 
                    WHERE id = ?
                ''',
                (gene_id.upper(),)
            )]
        else:
            cur = self._db.cursor()
            al_list = cur.executemany('''
                SELECT alias,id FROM aliases WHERE id = ?
                ''',
                [(id,) for id in gene_id]
            )
            als = dict()
            for al,id in al_list:
                if id in als:
                    als[id].append(al)
                else:
                    als[id] = [al]
            return als

    def remove_aliases(self):
        self._db.cursor().execute('DELETE FROM aliases;')

    def has_annotations(self):
        cur = self._db.cursor()
        cur.execute('SELECT count(*) FROM func;')
        return (int(cur.fetchone()[0]) > 0)
    
    def get_annotations(self,item):
        # Build the query from all the genes provided
        if isinstance(item,(set,list)):
            ls = "{}".format("','".join([str(x) for x in item]))
            single = False
        else:
            ls = item
            single = True
        query = "SELECT * FROM func WHERE id IN ('{}');".format(ls)
        
        # Run the query and turn the result into a list of tuples
        cur = self._db.cursor()
        cur.execute(query)
        annotes = cur.fetchall()
        
        # If a list of genes was passed in, return a dictionary of lists
        if not single:
            res = {}
            for id,desc in annotes:
                if id in res:
                    res[id].append(desc)
                else:
                    res[id] = [desc]
        
        # Otherwise just return the list annotations
        else:
            res = []
            for id,desc in annotes:
                res.append(desc)
        return res
    
    def export_annotations(self, filename=None, sep="\t"):
        '''
            Make a table of all functional annotations.
        '''
        # Find the default filename
        if filename == None:
            filename = self.name + '_func.tsv'
        
        # Pull them all from sqlite
        cur = self._db.cursor()
        cur.execute("SELECT * FROM func;")
        
        # Used pandas to save it
        df = pd.DataFrame(cur.fetchall(),columns=['gene','desc']).set_index('gene')
        df.to_csv(filename,sep=sep)
    
    def add_annotations(self, filename, sep="\t", gene_col=0, skip_cols=None):
        ''' 
            Imports Annotation relationships from a csv file. By default will
            assume gene names are first column

            Parameters
            ----------
            filename : str 
                The file containing the annotations
            sep : str (default: \\t)
                The delimiter for the columns in the annotation file
            gene_col : int (default: 0)
                The index of the column containing the gene IDs
            skip_cols : default:None
                Optional names of columns to drop before adding 
                annotations

            Returns
            -------
            None if successful

        '''
        # import from file, assume right now that in correct order
        tbl = pd.read_table(filename,sep=sep,dtype=object)
        idx_name = tbl.columns[gene_col]
        tbl[idx_name] = tbl[idx_name].str.upper()
        # Set thie index to be the specified gene column
        tbl.set_index(idx_name,inplace=True)        
        
        # Drop columns if we need to
        if skip_cols is not None:
            # removing certain columns
            tbl.drop(tbl.columns[skip_cols],axis=1,inplace=True)
        
        # Get rid of any genes not in the refence genome
        cur = self._db.cursor()
        cur.execute('SELECT id FROM loci;')
        rm = set(tbl.index.values) - set([id[0] for id in cur.fetchall()])
        tbl.drop(rm,axis=0,inplace=True)
        del rm, cur
        
        # One Annotation per row, drop the nulls and duplicates
        tbl = tbl.reset_index()
        tbl = pd.melt(tbl,id_vars=idx_name,var_name='col',value_name='desc')
        tbl.drop('col',axis=1,inplace=True)
        tbl.dropna(axis=0,inplace=True)
        tbl.drop_duplicates(inplace=True)
        
        # Run the transaction to throw them in there
        cur = self._db.cursor()
        try:
            cur.execute('BEGIN TRANSACTION')
            cur.executemany(
                'INSERT INTO func VALUES (?,?)'
                ,tbl.itertuples(index=False))
            cur.execute('END TRANSACTION')
        
        except Exception as e:
            self.log.info("import failed: {}",e)
            cur.execute('ROLLBACK')
        
        # Make sure the indices are built
        self._build_indices()

    def remove_annotations(self):
        self._db.cursor().execute('DELETE FROM func;')

    @classmethod
    def from_DataFrame(cls,df,name,description,build,organism,
            chrom_col='chrom',start_col='start',stop_col='stop',
            id_col='ID'):
        '''
            Imports a RefGen object from a CSV.
        '''
        self = cls.create(name,description,type='RefGen')
        self._global('build',build)
        self._global('organism','organism')
        genes = list()
        for i,row in df.iterrows():
            genes.append(
                Gene(
                    row[chrom_col],int(row[start_col]),int(row[stop_col]),
                    id=row[id_col],build=build,organism=organism
                ).update(dict(row.items()))
            )
        self.add_gene(genes)
        self._build_indices()
        return self


    def copy(self,name,description):
        '''
            Creates a copy of a refgen with a new name and description.

            Parameters
            ----------
                name : str
                    Name of the copy refgen
                description : str
                    Short description of the reference genome
            Returns
            -------
                co.RefGen object containing the same genes and 
                chromosomems as the original.
        '''
        copy = self.create(name,description,'RefGen')
        copy._global('build',self.build)
        copy._global('organism',self.organism)
        # Should have the same chromosomems
        for chrom in self.iter_chromosomes():
            copy.add_chromosome(chrom)
        # Should have the same gene list
        copy.add_gene(self.iter_genes(),refgen=self)
        copy._build_indices()
        return copy

    @classmethod
    def filtered_refgen(cls, name, description, refgen, gene_list):
        '''
            Copies from a previous instance of refgen, making sure
            each gene is within gene list
        '''
        self = cls.create(name,description,'RefGen')
        self._global('build',refgen.build)
        self._global('organism',refgen.organism)
        # Should have the same chromosomes
        for chrom in refgen.iter_chromosomes():
            self.add_chromosome(chrom)
        # Add the genes from genelist
        self.add_gene(gene_list,refgen=refgen)
        self._build_indices()
        return self


