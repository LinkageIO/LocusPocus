import matplotlib.pylab as plt
    
def plot_loci(self,loci,filename,flank_limit=2): # pragma: no cover
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
        cax.hold(True)
        # place marker for start window
        cax.scatter(hoffset,voffset,marker='>')
        # place marker for start snp
        cax.scatter(hoffset+locus.window,voffset,marker='.',color='blue')
        # place marker for stop snp
        cax.scatter(hoffset+locus.window+len(locus),voffset,marker='.',color='blue')
        # place marker for stop snp
        cax.scatter(hoffset+locus.window+len(locus)+locus.window,voffset,marker='<')

        # place markers for sub snps
        #for subsnp in locus.sub_loci:
        #    cax.scatter(
        #        hoffset + subsnp.start - locus.start + locus.window,
        #        voffset,
        #        marker='.',
        #        color='blue'
        #    )

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

