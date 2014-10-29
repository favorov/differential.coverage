#Differential.Coverage library
#A. Favorov, E. Fertig, S. Wheelan 2014
#annotation utilities

#'nucl.chromosomes.hg19
#'
#'Returns the choromosome lengths for hg19
#'
#'Returns the choromosome lengths for hg19 as a \code{seqinfo} object.
#'@param chrM boolean: if FALSE (default), returns info only for chr 1..22, X and Y. If TRUE, 1..22, X, Y and M
#'@return \code{seqinfo} object with the chromosome desriptions
#'@seealso seqinfo
nucl.chromosomes.hg19<-function(chrM=FALSE)
{
	chr.all<-seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)
	if (chrM)
		chr.all[names(chr.all)[1:25]]
	else
		chr.all[names(chr.all)[1:24]]
}


#'genes.by.TSS.overlap
#'
#'Generates list of genes that start inside a given set of intervals
#'
#'After the noodles (the set of intervals to search TSS in) are inflated by flanks, we look for all the TSS that start inside the (inflated)  intervals according to \code{TxDb.Hsapiens.UCSC.hg19.knownGene}. If the noodles has p.value and/or fdr metadata, we ascribe the data of the interval to the retrieved gene. If there a gene refers to a set of noodles, it has min ascribed. The ishyper data is also transferred to gene, if it is not contradictory.
#'@param noodles the list of intervals to look TSS in
#'@param flanks lenght to inflate the noddles by before the search
#'@return \code{GRanges} object that is the list of the genes we look for
genes.with.TSS.covered<-function(
	noodles, # GRanges with the noodles, if it has p.value, fdr and ishyper values, they will be mapped to genes
	flanks=0 #how far to shrink 
)
{
	expanded.noodles<-noodles
	#inflate DM noodles
	start(expanded.noodles)<-pmax(1,start(noodles)-flanks)
	end(expanded.noodles)<-pmin(end(noodles)+flanks,as.integer(seqlengths(noodles)[as.character(seqnames(noodles))]))
	#inflated

	#prepare gene TSS
	genelist<- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
	#initialise the list to subset later

	TSS<- genelist
	#this copy is for overlaps

	geneSymbols <- select(
		org.Hs.eg.db,
		keys=as.character(TSS$gene_id),
		columns=c('SYMBOL'),
		keytype='ENTREZID'
	)

	TSS$SYMBOL <- geneSymbols$SYMBOL

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
	

	noodle.TSS.Genes$SYMBOL<-TSS$SYMBOL[noodle.TSS.Gene.Indices]

	
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




#'gene.list.for.ovelapping.intervals
#'
#'Generates list of genes that start inside a given set of intervals
#'
#'After the noodles (the set of intervals to search TSS in) are inflated by flanks, we look for all the TSS that start inside the (inflated)  intervals according to \code{TxDb.Hsapiens.UCSC.hg19.knownGene}. 
#'If a noodle (interval) overlaps more that one TSS, we form a text list of the genes.
#'@inheritParams genes.with.TSS.covered
#'@return \code{GRanges} object, noodles argument with added genes for each interval  
genes.with.TSS.covered.by.interval<-function(
	noodles, # GRanges with the noodles, if it has p.value ans ishyper values, thay will be mapped to genes
	flanks=0 #how far to shrink 
)
{
	expanded.noodles<-noodles
	#inflate DM noodles
	start(expanded.noodles)<-pmax(1,start(noodles)-flanks)
	end(expanded.noodles)<-pmin(end(noodles)+flanks,as.integer(seqlengths(noodles)[as.character(seqnames(noodles))]))
	#prepare gene TSS
	genelist<- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

	TSS<- genelist

	geneSymbols <- select(
		org.Hs.eg.db,
		keys=as.character(TSS$gene_id),
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
	
	decorated.nooles<-noodles

	#overlapped.TSS

	overlapped.TSS<-tapply(TSS$SYMBOL[subjectHits(overla)],queryHits(overlapa),paste,collapse=', ')

	decorated.noodles<-cbind(decorated.noodles,'overlapped.TSS'=overlapped.TSS[as.character(1:length(decorated.noodles))])

	overlapped.pos<-tapply(as.character(start(TSS))[subjectHits(overla)],queryHits(overlapa),paste,collapse=', ')

	decorated.noodles<-cbind(decorated.noodles,'overlapped.pos'=overlapped.pos[as.character(1:length(decorated.noodles))])

	ovrl.dir<-tapply(as.character(strand(TSS))[subjectHits(overlapa)],queryHits(overlapa),paste,collapse=', ')

	decorated.noodles<-cbind(decorated.noodles,'ovrl.dir'=ovrl.dir[as.character(1:length(decorated.noodles))])

	message('mapped')

	decorated.noodles

}


