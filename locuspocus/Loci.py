#!/usr/bin/python3
import locuspocus.LocusDist as LocusDist

from collections import defaultdict

from .Locus import Gene,Locus
from .Exceptions import ZeroWindowError

import itertools
import collections
import random
import pandas as pd
import numpy as np
import scipy as sp
import math
import gzip
import re


class Loci(object):
    '''
        Just a bunch of Locuses.        
    '''
    def __init__(self):
        self.loci = pd.DataFrame(
            columns=[
                pd.Series([],name='chrom',dtype='object'),
                pd.Series([],name='start',dtype='int'),
                pd.Series([],name='end',dtype='int'),
                pd.Series([],name='id',dtype='object'),
                pd.Series([],name='window',dtype='int')
            ]
        )
        self.loci_attrs = pd.DataFrame()
        self.sub_loci = pd.DataFrame()

    def append(self,locus):
        if isinstance(locus,Locus):
            self.loci = self.loci.append(locus.as_dict(),ignore_index=True)
        else:
            # support adding iterables of loci
            loci = iter(locus)
            self.loci = self.loci.append(
                [locus.as_dict() for locus in loci],
                ignore_index=True
            )

    def __len__(self):
        '''
            Returns the number of genes in the dataset
        '''
        return len(self.loci)

    def __contains__(self,x):
        ''' flexible on what you pass into the 'in' function '''
        if isinstance(x,Locus):
            # find the locus by 
            if x.id in self.loci.id.values:
               return True
            else:
                return False
        elif isinstance(x,str):
            # Can be a string object
            if x in self.loci.id.values:
                return True
            else:
                return False
        else:
            raise TypeError('Cannot test for containment for {}'.format(obj))


#   def rand(self,**kwargs):
#       '''
#           Returns a random gene within the reference genome.
#           Also allows passing of keyword arguments to Locus
#           constructor method allowing for flexible generation.
#           See Locus.__init__ for more details.

#           Parameters
#           ----------
#           **kwargs : key,value pairs
#               Extra parameters passed onto the locus init method.

#           Returns
#           -------
#           A Gene object (camoco.Locus based)

#       '''
#       return self.Gene(*self.db.cursor().execute('''
#           SELECT chromosome,start,end,id from genes WHERE rowid = ?
#           ''',(random.randint(1,self.num_genes()),)).fetchone(),
#           **kwargs
#       )

#   def random_genes(self,n,**kwargs):
#       '''
#           Return random genes from the RefGen

#           Parameters
#           ----------
#           n : int

#           **kwargs : key,value pairs
#               Extra parameters passed onto the locus init method

#           Returns
#           -------
#           An iterable containing random genes

#       '''
#       rand_nums = np.random.randint(1,high=self.num_genes(),size=n)
#       gene_info = self.db.cursor().executemany(
#               "SELECT chromosome,start,end,id from genes WHERE rowid = ?",
#               [[int(rownum)] for rownum in rand_nums]
#       )
#       return set([Gene(chr,start,end=end,id=id,**kwargs) for \
#           (chr,start,end,id) in gene_info])


#   def iter_genes(self):
#       '''
#           Iterates over genes in RefGen.

#           Returns
#           -------
#           A generator containing genes
#       '''
#       for x in self.db.cursor().execute('''
#               SELECT chromosome,start,end,id FROM genes
#           '''):
#           yield self.Gene(*x,build=self.build,organism=self.organism)

#   def from_ids(self, gene_list, check_shape=False):
#       '''
#           Returns a list of gene object from an iterable of id strings
#           OR from a single gene id string.

#           Parameters
#           ----------
#           gene_list : str OR iterable of str
#               ID(s) of the genes you want to pull out
#           check_shape : bool (default: False)
#               Check if you get back the same number of ids you
#               pass in. If false (default), just give back what
#               you find, ignoring erronous ids.

#           Returns
#           -------
#           A list of locus objects if you pass in an iterable,
#           otherwise a single gene

#       '''

#       cur = self.db.cursor()
#       if isinstance(gene_list,str):
#       # Handle when we pass in a single id
#           gene_id = gene_list
#           if gene_id not in self:
#               result = cur.execute('SELECT id FROM aliases WHERE alias = ?', [gene_id]).fetchone()
#               if not result:
#                   raise ValueError('{} not in {}'.format(gene_id,self.name))
#               gene_id = result[0]
#           info = cur.execute('SELECT chromosome,start,end,id FROM genes WHERE id = ?', [gene_id]).fetchone()
#           return self.Gene(*info,build=self.build,organism=self.organism)

#       else:
#       # Handle when we pass an iterable of gene ids
#           bad_ids = []
#           gene_info = []
#           for id in gene_list:
#               gene_id = id
#               if gene_id not in self:
#                   result = cur.execute('SELECT id FROM aliases WHERE alias = ?', [gene_id]).fetchone()
#                   if not result:
#                       bad_ids.append(gene_id)
#                       continue
#                   gene_id = result[0]
#               gene_info.append(cur.execute('SELECT chromosome,start,end,id FROM genes WHERE id = ?', [gene_id]).fetchone())

#           genes = [self.Gene(*x,build=self.build,organism=self.organism) \
#                   for x in gene_info]

#           if check_shape and len(genes) != len(gene_list):
#               err_msg = '\nThese input ids do not have genes in reference:'
#               for x in range(len(bad_ids)):
#                   if x % 5 == 0:
#                       err_msg += '\n'
#                   err_msg += bad_ids[x] + '\t'
#               raise ValueError(err_msg)
#           return genes

#   def __getitem__(self,item):
#       '''
#           A convenience method to extract loci from the reference geneome.
#       '''
#       if isinstance(item,str):
#           # Handle when we pass in a single id
#           gene_id = item.upper()
#           try:
#               return self.Gene(
#                   *self.db.cursor().execute('''
#                       SELECT chromosome,start,end,id FROM genes WHERE id = ?
#                       ''',(gene_id,)
#                   ).fetchone(),
#                   build=self.build,
#                   organism=self.organism
#               )
#           except TypeError as e:
#               raise ValueError('{} not in {}'.format(gene_id,self.name))
#       genes = [
#           self.Gene(*x,build=self.build,organism=self.organism) \
#           for x in self.db.cursor().execute('''
#               SELECT chromosome,start,end,id FROM genes WHERE id IN ('{}')
#           '''.format("','".join(map(str.upper,item))))
#       ]
#       return genes

#   def genes_within(self,loci,chain=True):
#       '''
#           Returns the genes that START within a locus 
#           start/end boundry.

#           Looks like: (y=yes,returned; n=no,not returned)

#           nnn  nnnnnnn  yyyyyy   yyyyy  yyyyyy yyyyyyy
#                   start                        end
#               -----x****************************x-------------

#       '''
#       if isinstance(loci,Locus):
#           return [
#               self.Gene(*x,build=self.build,organism=self.organism) \
#               for x in self.db.cursor().execute('''
#                   SELECT chromosome,start,end,id FROM genes
#                   WHERE chromosome = ?
#                   AND start >= ? AND start <= ?
#                   ''',
#                   (loci.chrom,loci.start,loci.end))
#           ]
#       else:
#           iterator = iter(loci)
#           genes = [self.genes_within(locus,chain=chain) for locus in iterator]
#           if chain:
#               genes = list(itertools.chain(*genes))
#           return genes

#   def upstream_genes(self,locus,gene_limit=1000,window_size=None):
#       '''
#           Find genes that START upstream of a locus. 
#           Genes are ordered so that the nearest genes are 
#           at the beginning of the list.

#           Return Genes that overlap with the upstream window,
#           This includes partially overlapping genes, but NOT
#           genes that are returned by the genes_within method. 

#           Looks like: (y=yes,returned; n=no,not returned)

#           nnn  yyyyyyy   yyyyyy   yyyyy  yyyyyy nnnn nnnn nnnnnnnn
#                                            nnnn
#                                          start             end
#               -----------------------------x****************x--
#                  ^_________________________| Window (upstream)
#       '''
#       if locus.window == 0 and window_size is None:
#           raise ZeroWindowError(
#               'Asking for upstream genes for {}',
#               locus.id
#           )
#       if window_size is not None:
#           upstream = locus.start - window_size
#       else:
#           upstream = locus.upstream
#       return [
#           self.Gene(*x,build=self.build,organism=self.organism) \
#           for x in self.db.cursor().execute('''
#               SELECT chromosome,start,end,id FROM genes
#               WHERE chromosome = ?
#               AND start < ? -- Gene must start BEFORE locus
#               AND end >= ?  -- Gene must end AFTER locus window (upstream) 
#               ORDER BY start DESC
#               LIMIT ?
#           ''',(locus.chrom, locus.start, upstream, gene_limit)
#       )]

#   def downstream_genes(self,locus,gene_limit=1000,window_size=None):
#       '''
#           Returns genes downstream of a locus. Genes are ordered 
#           so that the nearest genes are at the beginning of the list.

#           Return Genes that overlap with the downstream window,
#           This includes partially overlapping genes, but NOT
#           genes that are returned by the genes_within method. 

#           Looks like: (y=yes,returned; n=no,not returned)

#           nnn  nnnnnnn   nnnnnn nnnn  yyyy  yyyyyy yyyy yyyyyy  nnnnn
#              start             end
#             ---x****************x--------------------------------
#                                 |_______________________^ Window (downstream)
#       '''
#       if locus.window == 0 and window_size is None:
#           raise ZeroWindowError(
#               'Asking for upstream genes for {} with no window size',
#               locus.id
#           )
#       if window_size is not None:
#           downstream = locus.end + window_size
#       else:
#           downstream = locus.downstream

#       return [
#           self.Gene(*x,build=self.build,organism=self.organism) \
#           for x in self.db.cursor().execute('''
#               SELECT chromosome,start,end,id FROM genes
#               WHERE chromosome = ?
#               AND start > ?
#               AND start <= ?
#               ORDER BY start ASC
#               LIMIT ?
#           ''',(locus.chrom, locus.end, downstream, gene_limit)
#       )]

#   def flanking_genes(self, loci, flank_limit=2,chain=True,window_size=None):
#       '''
#           Returns genes upstream and downstream from a locus
#           ** done NOT include genes within locus **
#       '''
#       if isinstance(loci,Locus):
#           # If we cant iterate, we have a single locus
#           locus = loci
#           if locus.window == 0 and window_size is None:
#               raise ZeroWindowError(
#                   'Asking for upstream genes for {} and no window size.',
#                   locus.id
#               )
#           upstream_gene_limit = int(flank_limit)
#           downstream_gene_limit = int(flank_limit)
#           up_genes = self.upstream_genes(
#               locus, gene_limit=upstream_gene_limit, window_size=window_size
#           )
#           down_genes = self.downstream_genes(
#               locus, gene_limit=downstream_gene_limit, window_size=window_size
#           )
#           if chain:
#               return list(itertools.chain(up_genes,down_genes))
#           return (up_genes,down_genes)
#       else:
#           iterator = iter(loci)
#           genes = [
#               self.flanking_genes(locus,flank_limit=flank_limit,window_size=window_size)\
#               for locus in iterator
#           ]
#           if chain:
#               genes = list(itertools.chain(*genes))
#           return genes

#   def candidate_genes(self, loci, flank_limit=2,
#       chain=True, window_size=None, include_parent_locus=False,
#       include_parent_attrs=False, include_num_intervening=False, 
#       include_rank_intervening=False, include_num_siblings=False,
#       attrs=None):
#       '''
#           Locus to Gene mapping.
#           Return Genes between locus start and stop, plus additional
#           flanking genes (up to flank_limit)

#           Parameters
#           ----------
#           loci : camoco.Locus (also handles an iterable containing Loci)
#               a camoco locus or iterable of loci
#           flank_limit : int (default : 2)
#               The total number of flanking genes **on each side**
#               considered a candidate surrounding a locus
#           chain : bool (default : true)
#               Calls itertools chain on results before returning
#           window_size : int (default: None)
#               Optional parameter used to extend or shorten a locus
#               window from which to choose candidates from. If None,
#               the function will resort to what is available in the
#               window attribute of the Locus.
#           include_parent_locus : bool (default: False)
#               Optional parameter which will update candidate genes
#               'attr' attribute with the id of the parent locus
#               which contains it.
#           include_parent_attrs : iterable (default: False)
#               Optional parameter to include attributes from the parent
#               locus. Parent locus attrs specified here will be included. 
#               If effective loci is > 1, the maximum value will be
#               included. E.g. - including the SNP effect size with 
#               candidate genes.
#           include_num_intervening : bool (default: False)
#               Optional argument which adds an attribute to each 
#               candidate genes containing the rank of each gene
#               as a function of distance away from the parent 
#               locus. (i.e. the closest candidate is 1 and the 
#               furthest candidate is n)
#           include_rank_intervening : bool (default: False)
#               Optional argument which adds the rank of each
#               candidatea as a funtion of distance from the parent
#               Locus. i.e. The closest gene is ranked 1 and the furthest
#               gene is ranked n.
#           include_num_siblings : bool (default: False)
#               Optional argument which adds an attribute to each
#               candidate gene containing the number of total 
#               candidates (siblings) identifies at the locus.
#           attrs : dict (default: None)
#               An optional dictionary which will be updated to each
#               candidate genes attr value.

#           Returns
#           -------
#           a list of candidate genes (or list of lists if chain is False)

#       '''
#       if isinstance(loci,Locus):
#           # If not an iterator, its a single locus
#           locus = loci
#           genes_within = self.genes_within(locus)
#           up_genes,down_genes = self.flanking_genes(
#               locus, flank_limit=flank_limit, chain=False,
#               window_size=window_size
#           )

#           # This always returns candidates together, if 
#           # you want specific up,within and down genes
#           # use the specific methods
#           genes = sorted(itertools.chain(up_genes,genes_within,down_genes))
#           #include the number of effective loci
#           if include_rank_intervening == True:
#               ranks = sp.stats.rankdata([locus-x for x in genes])
#           # Iterate through candidate genes and propagate the 
#           # parental info
#           for i,gene in enumerate(genes):
#               gene.update({'num_effective_loci':len(locus.sub_loci)})
#               # include parent locus id if thats specified
#               if include_parent_locus == True:
#                   gene.update({'parent_locus':locus.id})
#               if include_rank_intervening == True:
#                   gene.update({'intervening_rank':ranks[i]})
#               # update all the parent_attrs
#               if include_parent_attrs:
#                   for attr in include_parent_attrs:
#                       attr_name = 'parent_{}'.format(attr)
#                       gene.update({attr_name: locus[attr]})
#           if include_num_intervening == True:
#               num_down = 0
#               num_up = 0
#               # Sort the genes by their distance from the locus
#               genes = sorted(genes,key=lambda x: x-locus)
#               for gene in genes:
#                   if locus in gene:
#                       gene.update({'num_intervening':-1})
#                   elif gene > locus:
#                       gene.update({'num_intervening':num_down})
#                       num_down += 1
#                   elif gene < locus:
#                       gene.update({'num_intervening':num_up})
#                       num_up += 1

#       
#           if include_num_siblings == True:
#               for gene in genes:
#                   gene.update({'num_siblings':len(genes)})
#           if attrs is not None:
#               for gene in genes:
#                   gene.update(attrs)
#           return genes

#       else:
#           iterator = iter(sorted(loci))
#           genes = [
#               # This is becoming a pain in the ass
#               self.candidate_genes(
#                   locus, 
#                   flank_limit=flank_limit,
#                   chain=chain,
#                   window_size=window_size,
#                   include_parent_locus=include_parent_locus,
#                   include_parent_attrs=include_parent_attrs,
#                   include_num_intervening=include_num_intervening,
#                   include_rank_intervening=include_rank_intervening,
#                   include_num_siblings=include_num_siblings,
#                   attrs=attrs
#               ) for locus in iterator
#           ]
#           if chain:
#               genes = list(set(itertools.chain(*genes)))
#           return genes

#   def bootstrap_candidate_genes(self, loci, flank_limit=2,
#       chain=True, window_size=None, include_parent_locus=False):
#       '''
#           Returns candidate genes which are random, but conserves
#           total number of overall genes.

#           Parameters
#           ----------
#           loci : camoco.Locus (also handles an iterable containing Loci)
#               a camoco locus or iterable of loci
#           flank_limit : int (default : 2)
#               The total number of flanking genes **on each side**
#               considered a candidate surrounding a locus
#           chain : bool (default : true)
#               Calls itertools chain on results before returning,
#           include_parent_locus : bool (default: False)
#               Optional parameter which will update candidate genes
#               'attr' attribute with the id of the parent locus
#               which contains it.

#           Returns
#           -------
#           a list of candidate genes (or list of lists if chain is False)

#       '''
#       if isinstance(loci,Locus):
#           # We now have a single locus
#           locus = loci
#           # grab the actual candidate genes
#           num_candidates = len(
#               self.candidate_genes(
#                   locus, flank_limit=flank_limit,
#                   chain=True, window_size=window_size
#               )
#           )
#           if num_candidates == 0:
#               return []
#           # Snps a random genes from the genome
#           random_gene = self.rand()
#           # Snag the same number of candidates
#           random_candidates = self.upstream_genes(
#               random_gene, 
#               gene_limit=num_candidates,
#               window_size=10e100
#           )
#           if len(random_candidates) != num_candidates:
#               # somehow we hit the end of a chromosome
#               # or something, just recurse
#               random_candidates = self.bootstrap_candidate_genes(
#                   locus,flank_limit=flank_limit,chain=True
#               )
#           if include_parent_locus == True:
#               for gene in random_candidates:
#                   gene.update({'parent_locus':random_gene.id})
#           return random_candidates
#       else:
#           # Sort the loci so we can collapse down
#           locus_list = sorted(loci)
#           seen = set()
#           bootstraps = list()
#           target = self.candidate_genes(
#               locus_list,flank_limit=flank_limit,
#               chain=False,window_size=window_size
#           )
#           target_accumulator = 0
#           candidate_accumulator = 0
#           self.log('target: {}, loci: {}',len(target),len(locus_list))
#           for i,(locus,targ) in enumerate(zip(locus_list,target)):
#               # compare downstream of last locus to current locus
#               candidates = self.bootstrap_candidate_genes(
#                   locus, flank_limit=flank_limit, 
#                   chain=True, window_size=window_size,
#                   include_parent_locus=include_parent_locus
#               )
#               # If genes randomly overlap, resample
#               while len(seen.intersection(candidates)) > 0:
#                   candidates = self.bootstrap_candidate_genes(
#                       locus, flank_limit=flank_limit,
#                       window_size=window_size, chain=True,
#                       include_parent_locus=include_parent_locus
#                   )
#               # Add all new bootstrapped genes to the seen list
#               seen |= set(candidates)
#               bootstraps.append(candidates)
#           if chain:
#               bootstraps = list(seen)
#           self.log("Found {} bootstraps",len(bootstraps))
#           return bootstraps


#   def pairwise_distance(self, gene_list=None):
#       '''
#           returns a vector containing the pairwise distances between genes
#           in gene_list in vector form. See np.squareform for matrix
#           conversion.
#       '''
#       if gene_list is None:
#           gene_list = list(self.iter_genes())
#       query = '''
#               SELECT genes.id, chrom.rowid, start, end FROM genes
#               LEFT JOIN chromosomes chrom ON genes.chromosome = chrom.id
#               WHERE genes.id in ("{}")
#               ORDER BY genes.id
#       '''.format('","'.join([g.id for g in gene_list]))
#       # extract chromosome row ids and gene start positions for each gene
#       positions = pd.DataFrame(
#           # Grab the chromosomes rowid because its numeric
#           self.db.cursor().execute(query).fetchall(),
#           columns=['gene','chrom','start','end']
#       ).sort('gene')
#       # chromosome needs to be floats
#       positions.chrom = positions.chrom.astype('float')
#       # Do a couple of checks
#       assert len(positions) == len(gene_list), \
#           'Some genes in dataset not if RefGen'
#       assert all(positions.gene == [g.id for g in gene_list]), \
#           'Genes are not in the correct order!'
#       distances = LocusDist.gene_distances(
#           positions.chrom.values,
#           positions.start.values,
#           positions.end.values
#       )
#       return distances

#   def summary(self):
#       print ("\n".join([
#           'Reference Genome: {} - {} - {}',
#           '{} genes']).format(
#               self.organism,self.build,
#               self.name,self.num_genes()
#           )
#       )

#   def __contains__(self,obj):
#       ''' flexible on what you pass into the 'in' function '''
#       if isinstance(obj,Locus):
#           # you can pass in a gene object (this expression
#           # should ALWAYS be true if you
#           # created gene object from this RefGen)
#           if self.db.cursor().execute(
#               '''SELECT COUNT(*) FROM genes WHERE id = ?''',
#               (obj.id.upper(),)).fetchone()[0] == 1:
#               return True
#           else:
#               return False
#       elif isinstance(obj,str):
#           # Can be a string object
#           if self.db.cursor().execute('''
#               SELECT COUNT(*) FROM genes WHERE id = ?''',
#               (str(obj).upper(),)).fetchone()[0] == 1:
#               return True
#           else:
#               return False
#       else:
#           raise TypeError('Cannot test for containment for {}'.format(obj))

#   def _build_indices(self):
#       self.log('Building Indices')
#       cur = self.db.cursor()
#       cur.execute('''
#           CREATE INDEX IF NOT EXISTS genepos ON genes (chromosome,start);
#           CREATE INDEX IF NOT EXISTS geneid ON genes (id);
#           CREATE INDEX IF NOT EXISTS geneattr ON gene_attrs (id);
#       ''')

#   def add_chromosome(self,chrom):
#       ''' adds a chromosome object to the class '''
#       self.db.cursor().execute('''
#           INSERT OR REPLACE INTO chromosomes VALUES (?,?)
#       ''',(chrom.id,chrom.length))

#   def add_aliases(self, alias_file, id_col=0, alias_col=1, headers=True):
#       ''' 
#           Add alias map to the RefGen 
#           
#           Parameters
#           ----------
#           alias_file : string (path)
#               The path to the alias file
#           id_col : int (default: 0)
#               The column containing the gene identifier
#           alias_col : int (default: 1)
#               The columns containing the alias
#           header : bool (default: True)
#               A switch stating if there is a header row to ignore or not
#       '''
#       with rawFile(alias_file) as IN:
#           if headers:
#               garb = IN.readline()
#   
#           aliases = []
#           self.log('Importing aliases from: {}',alias_file)
#           for line in IN:
#               row = re.split(',|\t',line)
#               if row[id_col].strip() in self:
#                   aliases.append(
#                       (row[alias_col],row[id_col].strip())
#                   )
#       cur = self.db.cursor()
#       self.log('Saving them in the alias table.')
#       cur.execute('BEGIN TRANSACTION')
#       cur.executemany(
#           'INSERT OR REPLACE INTO aliases VALUES (?,?)',
#           aliases
#       )
#       cur.execute('END TRANSACTION')

#   def num_aliases(self):
#       '''
#           Returns the number of aliases currently in the database
#       '''
#       return self.db.cursor() \
#           .execute('SELECT COUNT(*) FROM aliases') \
#           .fetchone()[0]

#   def aliases(self, gene_id):
#       if isinstance(gene_id,str):
#           return [alias[0] for alias in self.db.cursor().execute('''
#                   SELECT alias FROM aliases 
#                   WHERE id = ?
#               ''',
#               (gene_id.upper(),)
#           )]
#       else:
#           cur = self.db.cursor()
#           al_list = cur.executemany('''
#               SELECT alias,id FROM aliases WHERE id = ?
#               ''',
#               [(id,) for id in gene_id]
#           )
#           als = dict()
#           for al,id in al_list:
#               if id in als:
#                   als[id].append(al)
#               else:
#                   als[id] = [al]
#           return als

#   @classmethod
#   def from_gff(cls,filename,name,description,build,organism,
#                chrom_feature='chromosome',gene_feature='gene',
#                ID_attr='ID',attr_split='='):
#       '''
#           Imports RefGen object from a gff (General Feature Format) file.
#           See more about the format here:
#           http://www.ensembl.org/info/website/upload/gff.html

#           Parameters
#           ----------

#           filename : str
#               The path to the GFF file.
#           name : str
#               The name if the RefGen object to be stored in the core
#               camoco database.
#           description : str
#               A short description of the RefGen for future reference
#           build : str
#               A string designating the genome build, used for comparison
#               operations, genes may share IDS but are different across build.
#           organism : str
#               A short string describing the organims this RefGen is coming
#               from. Again, is used in comparing equality among genes which
#               may have the same id or name.
#           chrom_feature : str (default: chromosome)
#               The name of the feature (in column 3) that designates a
#               a chromosome.
#           gene_feature : str (default: gene)
#               The name of the feature (in column 2) that designates a 
#               gene. These features will be the main object that the RefGen
#               encompasses. 
#           ID_attr : str (default: ID)
#               The key in the attribute column which designates the ID or 
#               name of the feature.
#           attr_split : str (default: '=')
#               The delimiter for keys and values in the attribute column

#       '''
#       self = cls.create(name,description,type='RefGen')
#       self._global('build',build)
#       self._global('organism',organism)
#       genes = list()
#       chroms = dict()
#       if filename.endswith('.gz'):
#           IN = gzip.open(filename,'rt')
#       else:
#           IN = open(filename,'r')
#       for line in IN:
#           #skip comment lines
#           if line.startswith('#'):
#               continue
#           (chrom,source,feature,start,
#            end,score,strand,frame,attributes) = line.strip().split('\t')
#           attributes = dict([(field.strip().split(attr_split)) \
#               for field in attributes.strip(';').split(';')])
#           if feature == gene_feature:
#               genes.append(
#                   Gene(
#                       chrom,int(start),int(end),
#                       attributes[ID_attr].upper().strip('"'),strand=strand,
#                       build=build,organism=organism,**attributes
#                   ).update(attributes)
#               )
#               # Keep track of seen chromosomes
#               if chrom not in chroms:
#                   chroms[chrom] = end
#               else:
#                   if end > chroms[chrom]:
#                       chroms[chrom] = end
#       IN.close()
#       self.add_gene(genes)
#       self._build_indices()
#       return self

#   @classmethod
#   def filtered_refgen(cls,name,description,refgen,gene_list):
#       '''
#           Copies from a previous instance of refgen, making sure
#           each gene is within gene list
#       '''
#       self = cls.create(name,description,'RefGen')
#       self._global('build',refgen.build)
#       self._global('organism',refgen.organism)
#       # Should have the same chromosomes
#       for chrom in refgen.iter_chromosomes():
#           self.add_chromosome(chrom)
#       # Add the genes from genelist
#       self.add_gene(gene_list,refgen=refgen)
#       self._build_indices()
#       return self

#   def _create_tables(self):
#       cur = self.db.cursor()
#       cur.execute('''
#           CREATE TABLE IF NOT EXISTS chromosomes (
#               id TEXT NOT NULL UNIQUE,
#               length INTEGER NOT NULL
#           );
#           CREATE TABLE IF NOT EXISTS genes (
#               id TEXT NOT NULL UNIQUE,
#               chromosome TEXT NOT NULL,
#               start INTEGER,
#               end INTEGER
#           );
#           CREATE TABLE IF NOT EXISTS gene_attrs (
#               id TEXT NOT NULL,
#               key TEXT,
#               val TEXT
#           );
#           CREATE TABLE IF NOT EXISTS aliases (
#               alias TEXT UNIQUE,
#               id TEXT
#           );''');
