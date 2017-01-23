#differential.coverage library
#A. Favorov, E. Fertig, D.Gaykalova, J. Califano, S. Wheelan 2014-2016
#annotation utilities

#get the name of the TxDb object by the genome name
.knownGenes.by.genome.id<-function(genome.id)
{
	if (genome.id=='hg38')
		return('TxDb.Hsapiens.UCSC.hg38.knownGene')
	if (genome.id=='hg19')
		return('TxDb.Hsapiens.UCSC.hg19.knownGene')
	if(genome.id=='hg18')
		return('TxDb.Hsapiens.UCSC.hg18.knownGene')
	stop(paste0('I cannot make annotation for genome ',genome.id))
}


.getKnownGeneList<-function(genome.id)
{
	#prepare genes; we refere the TxDb object by name
	genelist<-unlist(
		genes(
			get(.knownGenes.by.genome.id(genome.id)),
			single.strand.genes.only=FALSE
		)
	)
	#we neeeded gene_id field for each gene
	#genelist$gene_id=names(genelist)
	#actually, names() is enough
	genelist
}

#'genes.with.TSS.covered
#'
#'Generates list of genes that start inside a given set of intervals
#'
#'After the noodles (the set of intervals to search TSS in) are inflated by flanks, we look for all the TSS that start inside the (inflated)  intervals according to \code{TxDb.Hsapiens.UCSC.hg19.knownGene} for \code{genome.id==hg19} (default), acording to \code{TxDb.Hsapiens.UCSC.hg18.knownGene} for \code{genome.id==hg18} and to \code{TxDb.Hsapiens.UCSC.hg38.knownGene} for \code{genome.id==hg38}. If the noodles has p.value and/or fdr metadata, we ascribe the data of the interval to the retrieved gene. If there a gene refers to a set of noodles, it has min ascribed. The ishyper data is also transferred to gene, if it is not contradictory.
#'@param noodles the \code{GRanges} list of intervals to look TSS in
#'@param flanks lenght to inflate the noddles by before the search; if >0, the seqlenght information is to be set in \code{noodles}
#'@param genome.id the character string with the id of genome we work with, the default is 'hg19', currntly, we work with hg19 or hg18
#'@return \code{GRanges} object that is the list of the genes we look for - the object is not co-indexed with \code{noodles} parameter
genes.with.TSS.covered<-function(
	noodles, # GRanges with the noodles, if it has p.value, fdr and ishyper values, they will be mapped to genes
	flanks=0, #how far to shrink
	genome.id='hg19' 
)
{
	knownGenes<-.knownGenes.by.genome.id(genome.id)
	expanded.noodles<-noodles
	#inflate DM noodles
	start(expanded.noodles)<-pmax(1,start(noodles)-flanks)
	end(expanded.noodles)<-pmin(end(noodles)+flanks,as.integer(seqlengths(noodles)[as.character(seqnames(noodles))]))
	#inflated

	#prepare genes; we refere the TxDb object by name
	TSS<-.getKnownGeneList(genome.id)
	#initialise the list to subset later

	#prepare TSS
	TSS<-genelist
	#this copy is for overlaps

	geneSymbols <- select(
		org.Hs.eg.db,
		keys=as.character(names(TSS)),
		columns=c('SYMBOL'),
		keytype='ENTREZID'
	)

	genelist$SYMBOL <- geneSymbols$SYMBOL

	tss.start<-ifelse(strand(TSS)=='+',start(TSS),end(TSS))

	start(TSS)<-tss.start
	end(TSS)<-tss.start
	#TSS prepared

	#make overlap
	overlaps<-findOverlaps(expanded.noodles,TSS)
	message('overlapped')

	noodle.TSS.Gene.Indices<-unique(subjectHits(overlaps))
	#indices of all the genes tha were hit by any noodle
	#maybe, some of them are hit with more than one

	noodle.TSS.Genes<-genelist[noodle.TSS.Gene.Indices]
	

	noodle.TSS.Genes$SYMBOL<-genelist$SYMBOL[noodle.TSS.Gene.Indices]

	
	if ('ishyper' %in% names(elementMetadata(noodles)))
	{
		ishyper<- tapply(noodles[queryHits(overlaps)]$ishyper,subjectHits(overlaps),
			function(li)
			{	
				vals<-unique(li)
				if (length(vals)==1)
					vals[1]
				else NA
			}
		)
		noodle.TSS.Genes$ishyper<-ishyper[as.character(noodle.TSS.Gene.Indices)]
		# we make the addressin because tapply return keys sorted, nothin to do with original order
		#and, we address by kyes
	}

	if ('p.value' %in% names(elementMetadata(noodles)))
	{
		p.value<-tapply(noodles[queryHits(overlaps)]$p.value,subjectHits(overlaps),min)
		# print(p.value)
		noodle.TSS.Genes$p.value<-p.value[as.character(noodle.TSS.Gene.Indices)] 
		# we make the addressin because tapply return keys sorted, nothin to do with original order
	}

	if ('fdr' %in% names(elementMetadata(noodles)))
	{
		fdr<-tapply(noodles[queryHits(overlaps)]$fdr,subjectHits(overlaps),min)
		# print(p.value)
		noodle.TSS.Genes$fdr<-fdr[as.character(noodle.TSS.Gene.Indices)] 
		# we make the addressin because tapply return keys sorted, nothin to do with original order
	}
	
	if ('FDR' %in% names(elementMetadata(noodles)))
	{
		FDR<-tapply(noodles[queryHits(overlaps)]$FDR,subjectHits(overlaps),min)
		# print(p.value)
		noodle.TSS.Genes$FDR<-FDR[as.character(noodle.TSS.Gene.Indices)] 
		# we make the addressin because tapply return keys sorted, nothin to do with original order
	}

	message('mapped')
	noodle.TSS.Genes
}

#'genes.intersected
#'
#'Generates list of genes that intersect a given set of intervals
#'
#'After the noodles (the set of intervals to search TSS in) are inflated by flanks, we look for all the genes that intersect the (inflated)  intervals according to \code{TxDb.Hsapiens.UCSC.hg19.knownGene} for \code{genome.id==hg19} (default), acording to \code{TxDb.Hsapiens.UCSC.hg18.knownGene} for \code{genome.id==hg18} and to \code{TxDb.Hsapiens.UCSC.hg38.knownGene} for \code{genome.id==hg38}. If the noodles has p.value and/or fdr metadata, we ascribe the data of the interval to the retrieved gene. If there a gene refers to a set of noodles, it has min ascribed. The ishyper data is also transferred to gene, if it is not contradictory.
#'@inheritParams genes.with.TSS.covered
#'@return \code{GRanges} object that is the list of the genes we look for - the object is not co-indexed with \code{noodles} parameter
genes.intersected<-function(
	noodles, # GRanges with the noodles, if it has p.value, fdr and ishyper values, they will be mapped to genes
	flanks=0, #how far to shrink
	genome.id='hg19' 
)
{
	expanded.noodles<-noodles
	#inflate DM noodles
	start(expanded.noodles)<-pmax(1,start(noodles)-flanks)
	end(expanded.noodles)<-pmin(end(noodles)+flanks,as.integer(seqlengths(noodles)[as.character(seqnames(noodles))]))
	#inflated

	#prepare genes; we refere the TxDb object by name

	genelist<-.getKnownGeneList(genome.id)
	#initialise the list to subset later

	geneSymbols <- select(
		org.Hs.eg.db,
		keys=as.character(names(genelist)),
		columns=c('SYMBOL'),
		keytype='ENTREZID'
	)

	genelist$SYMBOL <- geneSymbols$SYMBOL


	#make overlap
	overlaps<-findOverlaps(expanded.noodles,genelist)
	message('overlapped')

	noodle.Gene.Indices<-unique(subjectHits(overlaps))
	#indices of all the genes tha were hit by any noodle
	#maybe, some of them are hit with more than one

	noodle.Genes<-genelist[noodle.Gene.Indices]
	

	noodle.Genes$SYMBOL<-genelist$SYMBOL[noodle.Gene.Indices]

	
	if ('ishyper' %in% names(elementMetadata(noodles)))
	{
		ishyper<- tapply(noodles[queryHits(overlaps)]$ishyper,subjectHits(overlaps),
			function(li)
			{	
				vals<-unique(li)
				if (length(vals)==1)
					vals[1]
				else NA
			}
		)
		noodle.Genes$ishyper<-ishyper[as.character(noodle.Gene.Indices)]
		# we make the addressin because tapply return keys sorted, nothin to do with original order
		#and, we address by kyes
	}

	if ('p.value' %in% names(elementMetadata(noodles)))
	{
		p.value<-tapply(noodles[queryHits(overlaps)]$p.value,subjectHits(overlaps),min)
		# print(p.value)
		noodle.Genes$p.value<-p.value[as.character(noodle.Gene.Indices)] 
		# we make the addressin because tapply return keys sorted, nothin to do with original order
	}

	if ('fdr' %in% names(elementMetadata(noodles)))
	{
		fdr<-tapply(noodles[queryHits(overlaps)]$fdr,subjectHits(overlaps),min)
		# print(p.value)
		noodle.Genes$fdr<-fdr[as.character(noodle.Gene.Indices)] 
		# we make the addressin because tapply return keys sorted, nothin to do with original order
	}
	
	if ('FDR' %in% names(elementMetadata(noodles)))
	{
		FDR<-tapply(noodles[queryHits(overlaps)]$FDR,subjectHits(overlaps),min)
		# print(p.value)
		noodle.Genes$FDR<-FDR[as.character(noodle.Gene.Indices)] 
		# we make the addressin because tapply return keys sorted, nothin to do with original order
	}

	message('mapped')
	noodle.Genes
}


#'genes.with.TSS.covered.by.interval
#'
#'Generates a list of genes (possibly, empty) that start inside each interval
#'
#'After the noodles (the set of intervals to search TSS in) are inflated by flanks, we look for all the TSS that start inside the (inflated)  intervals according to \code{TxDb} object we use (TxDb.Hsapiens.UCSC.hg19.knownGene for genome.id=='hg19', TxDb.Hsapiens.UCSC.hg18.knownGene for hg18 or TxDb.Hsapiens.UCSC.hg38.knownGene for hg38). 
#'If a noodle (interval) overlaps more that one TSS, we form a text list of the genes.
#'@inheritParams genes.with.TSS.covered
#'@return \code{GRanges} object, noodles argument with added TSS-overlapped genes for each interval  
genes.with.TSS.covered.by.interval<-function(
	noodles, # GRanges with the noodles, if it has p.value ans ishyper values, thay will be mapped to genes
	flanks=0, #how far to shrink 
	genome.id='hg19' 
)
{
	knownGenes<-.knownGenes.by.genome.id(genome.id)
	expanded.noodles<-noodles
	#inflate DM noodles
	start(expanded.noodles)<-pmax(1,start(noodles)-flanks)
	end(expanded.noodles)<-pmin(end(noodles)+flanks,as.integer(seqlengths(noodles)[as.character(seqnames(noodles))]))
	#inflated
	#prepare gene TSS; we refere the TxDb object by name

	TSS<-.getKnownGeneList(genome.id)

	geneSymbols <- select(
		org.Hs.eg.db,
		keys=as.character(names(TSS)),
		columns=c('SYMBOL'),
		keytype='ENTREZID'
	)

	TSS$SYMBOL <- geneSymbols$SYMBOL

	tss.start<-ifelse(strand(TSS)=='+',start(TSS),end(TSS))

	start(TSS)<-tss.start
	end(TSS)<-tss.start
	#IAMHERE
	#make overlap
	overlapa<-findOverlaps(expanded.noodles,TSS)
	message('overlapped')
	
	decorated.noodles<-noodles

	#overlapped.TSS

	overlapped.TSS<-tapply(TSS$SYMBOL[subjectHits(overlapa)],queryHits(overlapa),paste,collapse=', ')
	decorated.noodles$overlapped.TSS=overlapped.TSS[as.character(1:length(decorated.noodles))]
	#we need this addressing scheme ([as.character(1:length(decorated.noodles))]) because names(overlapped.TSS)
	#are the indices of noodles that have overlapped TSS. Those that have not are not represented in overlapped.TSS
	#and after our indexing they are NA, and it is exacltly what we want them to be

	overlapped.pos<-tapply(as.character(start(TSS))[subjectHits(overlapa)],queryHits(overlapa),paste,collapse=', ')
	decorated.noodles$overlapped.pos=overlapped.pos[as.character(1:length(decorated.noodles))]

	ovrl.dir<-tapply(as.character(strand(TSS))[subjectHits(overlapa)],queryHits(overlapa),paste,collapse=', ')
	decorated.noodles$ovrl.dir=ovrl.dir[as.character(1:length(decorated.noodles))]

	message('mapped')

	decorated.noodles

}


#'genes.intersected.by.interval
#'
#'Generates a list of genes (possibly, empty) that intersects with each interval
#'
#'After the noodles (the set of intervals to search intersections with) are inflated by flanks, we look for all the genes that overlap with the (inflated)  intervals according to \code{TxDb} object we use (TxDb.Hsapiens.UCSC.hg19.knownGene for genome.id=='hg19', TxDb.Hsapiens.UCSC.hg18.knownGene for hg18 or TxDb.Hsapiens.UCSC.hg38.knownGene for hg38). 
#'If a noodle (interval) overlaps more that one gene, we form a text list of the genes.
#'@inheritParams genes.with.TSS.covered
#'@return \code{GRanges} object, noodles argument with added gene-overlapped genes for each interval  
genes.intersected.by.interval<-function(
	noodles, # GRanges with the noodles, if it has p.value ans ishyper values, thay will be mapped to genes
	flanks=0, #how far to shrink 
	genome.id='hg19' 
)
{
	knownGenes<-.knownGenes.by.genome.id(genome.id)
	expanded.noodles<-noodles
	#inflate DM noodles
	start(expanded.noodles)<-pmax(1,start(noodles)-flanks)
	end(expanded.noodles)<-pmin(end(noodles)+flanks,as.integer(seqlengths(noodles)[as.character(seqnames(noodles))]))
	#inflated
	#prepare gene list; we refere the TxDb object by name
	genelist<-.getKnownGeneList(genome.id)

	geneSymbols <- select(
		org.Hs.eg.db,
		keys=as.character(names(genelist)),
		columns=c('SYMBOL'),
		keytype='ENTREZID'
	)

	genelist$SYMBOL <- geneSymbols$SYMBOL

	overlapa<-findOverlaps(expanded.noodles,genelist)
	message('overlapped')
	
	decorated.noodles<-noodles

	#overlapped.genes

	overlapped.genes<-tapply(genelist$SYMBOL[subjectHits(overlapa)],queryHits(overlapa),paste,collapse=', ')
	decorated.noodles$overlapped.genes=overlapped.genes[as.character(1:length(decorated.noodles))]
	#we need this addressing scheme ([as.character(1:length(decorated.noodles))]) because names(overlapped.TSS)
	#are the indices of noodles that have overlapped TSS. Those that have not are not represented in overlapped.TSS
	#and after our indexing they are NA, and it is exacltly what we want them to be

	overlapped.pos<-tapply(as.character(start(genelist))[subjectHits(overlapa)],queryHits(overlapa),paste,collapse=', ')
	decorated.noodles$overlapped.pos=overlapped.pos[as.character(1:length(decorated.noodles))]

	ovrl.dir<-tapply(as.character(strand(genelist))[subjectHits(overlapa)],queryHits(overlapa),paste,collapse=', ')
	decorated.noodles$ovrl.dir=ovrl.dir[as.character(1:length(decorated.noodles))]

	message('mapped')

	decorated.noodles

}


#'closest.gene.by.interval
#'
#'Find the closest gene for each of the given set of intervals
#'
#'@inheritParams genes.with.TSS.covered
#'@return \code{GRanges} object, noodles argument with added closest gene info for each interval  
closest.gene.by.interval<-function(
	noodles, # GRanges with the noodles, if it has p.value ans ishyper values, thay will be mapped to genes
	genome.id='hg19' 
)
{
	knownGenes<-.knownGenes.by.genome.id(genome.id)
	message('closest')

	#prepare gene TSS; we refere the TxDb object by name
	genelist<-.getKnownGeneList(genome.id)

	geneSymbols <- select(
		org.Hs.eg.db,
		keys=as.character(names(genelist)),
		columns=c('SYMBOL'),
		keytype='ENTREZID'
	)

	genelist$SYMBOL <- geneSymbols$SYMBOL

	#now, TSS contains the gene start
  
	near.gene<-nearest(noodles,genelist)

	#some nearest contatain NA
	#we just remove them from result
	#mapping to no-NA indices
	is.a.near.gene <- !is.na(near.gene)

	noodles.nna<-noodles[is.a.near.gene]
	#decorated.noodles.nna now are all the noodles with non-NA nearest

	near.gene<-near.gene[is.a.near.gene]
	#near.gene now are all the non-NA values

	dist.gene<-distance(noodles.nna,genelist[near.gene])

	noodles.decoration<-DataFrame(
		closest.gene=character(length(noodles)),
		start=integer(length(noodles)),
		end=integer(length(noodles)),
		strand=character(length(noodles)),
		dist=integer(length(noodles))
	)

	noodles.decoration$closest.gene[is.a.near.gene]<-genelist$SYMBOL[near.gene]
	noodles.decoration$start[is.a.near.gene]<-start(genelist)[near.gene]
	noodles.decoration$end[is.a.near.gene]<-end(genelist)[near.gene]
	noodles.decoration$strand[is.a.near.gene]<-
		as.character(strand(genelist)[near.gene])
	noodles.decoration$dist[is.a.near.gene]<-dist.gene
	
	noodles.decoration[!is.a.near.gene,]=NA

	decorated.noodles<-noodles

	mcols(decorated.noodles)<-noodles.decoration

	message('mapped')
	
	decorated.noodles
}


#'closest.gene.start.by.interval
#'
#'Find the closest gene start for each of the given set of intervals
#'
#'@inheritParams genes.with.TSS.covered
#'@return \code{GRanges} object, noodles argument with added closest gene info for each interval  
closest.gene.start.by.interval<-function(
	noodles, # GRanges with the noodles, if it has p.value ans ishyper values, thay will be mapped to genes
	genome.id='hg19' 
)
{
	knownGenes<-.knownGenes.by.genome.id(genome.id)
	message('closest')

	#prepare gene TSS; we refere the TxDb object by name
	TSS<-.getKnownGeneList(genome.id)

	geneSymbols <- select(
		org.Hs.eg.db,
		keys=as.character(names(TSS)),
		columns=c('SYMBOL'),
		keytype='ENTREZID'
	)

	TSS$SYMBOL <- geneSymbols$SYMBOL

	tss.start<-ifelse(strand(TSS)=='+',start(TSS),end(TSS))

	start(TSS)<-tss.start
	end(TSS)<-tss.start

	#now, TSS contains the gene start
  
	near.TSS<-nearest(noodles,TSS)

	#some nearest contatain NA
	#we just remove them from result
	#mapping to no-NA indices
	is.a.near.TSS <- !is.na(near.TSS)

	noodles.nna<-noodles[is.a.near.TSS]
	#decorated.noodles.nna now are all the noodles with non-NA nearest

	near.TSS<-near.TSS[is.a.near.TSS]
	#near.TSS now are all the non-NA values

	dist.TSS<-distance(noodles.nna,TSS[near.TSS])


	dist.TSS<-ifelse(strand(TSS)[near.TSS]=='+',
			ifelse(start(noodles.nna)>start(TSS)[near.TSS],dist.TSS,-dist.TSS),
			ifelse(start(noodles.nna)>start(TSS)[near.TSS],-dist.TSS,dist.TSS)
	)

	noodles.decoration<-DataFrame(
		closest.TSS=character(length(noodles)),
		pos=integer(length(noodles)),
		dir=character(length(noodles)),
		dist=integer(length(noodles))
	)

	noodles.decoration$closest.TSS[is.a.near.TSS]<-TSS$SYMBOL[near.TSS]
	noodles.decoration$pos[is.a.near.TSS]<-start(TSS)[near.TSS]
	noodles.decoration$dir[is.a.near.TSS]<-as.character(strand(TSS)[near.TSS])
	noodles.decoration$dist[is.a.near.TSS]<-dist.TSS
	
	noodles.decoration[!is.a.near.TSS,]=NA

	decorated.noodles<-noodles

	mcols(decorated.noodles)<-noodles.decoration

	message('mapped')
	
	decorated.noodles
}
