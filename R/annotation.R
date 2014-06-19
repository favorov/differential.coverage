#Differential.Coverage library
#A. Favorov, E. Fertig, S. Wheelan 2014
#annotation utilities

nucl.chromosomes.hg19<-function()
{
	chr.all<-seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)
	chr.all[names(chr.all)[1:24]]
	#only chrNN, chrX and chrY.
	#chrM is 25
}


#we get a list of intervals with p-values and ishyper
#we inflate them
#we make an overall list of genes
gene.list.by.overlap<-function(
	noodles, # GRanges with the noodles, if it has p.value, fdr and ishyper values, thay will be mapped to genes
	flanks=0 #how far to shrink 
)
{
	expanded.noodles<-noodles
	#inflate DM noodles
	start(expanded.noodles)<-pmax(0,start(noodles)-flanks)
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

	#make overlap
	overlaps<-findOverlaps(expanded.noodles,TSS)
	message('overlapped')

	DM.Gene.Indices<-unique(subjectHits(overlaps))

	DM.Genes<-genelist[DM.Gene.Indices]

	DM.Genes$SYMBOL<-TSS$SYMBOL[DM.Gene.Indices]

	#print(DM.Gene.Indices)

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
		#print(ishyper)
		DM.Genes$ishyper<-ishyper[as.character(DM.Gene.Indices)]
		# we make the addressin because tapply return keys sorted, nothin to do with original order
	}

	if ('p.value' %in% names(elementMetadata(noodles)))
	{
		p.value<- tapply(noodles[queryHits(overlaps)]$p.value,subjectHits(overlaps),min)
		#print(p.value)
		DM.Genes$p.value<-p.value[as.character(DM.Gene.Indices)] 
		# we make the addressin because tapply return keys sorted, nothin to do with original order
	}

	message('mapped')

	DM.Genes
}



#we get a list of intervals 
#we inflate them
#we make a list of genes with overlapping TSS for each of them

gene.list.by.ovelapping.intervals<-function(
	noodles, # GRanges with the noodles, if it has p.value ans ishyper values, thay will be mapped to genes
	flanks=0 #how far to shrink 
)
{
	expanded.noodles<-noodles
	#inflate DM noodles
	start(expanded.noodles)<-pmax(0,start(noodles)-flanks)
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

	#make overlap
	overlaps<-findOverlaps(expanded.noodles,TSS)
	message('overlapped')
#I am here
	DM.CpGIs.GRanges<-as(DM.CpGIs.Ranges,'GRanges') # the same tester again

	end(DM.CpGIs.GRanges)<-pmin(end(DM.CpGIs.GRanges)+flanks,seqlengths(TSS)[as.character(seqnames(DM.CpGIs.GRanges))])

	start(DM.CpGIs.GRanges)<-pmax(start(DM.CpGIs.GRanges)-flanks,1)

	overla<-findOverlaps(DM.CpGIs.GRanges,TSS)

	#overlapped.TSS

	overlapped.TSS<-tapply(TSS$SYMBOL[subjectHits(overla)],queryHits(overla),paste,collapse=', ')

	DM.CpGIs.stat<-cbind(DM.CpGIs.stat,'overlapped.TSS'=overlapped.TSS[as.character(1:length(DM.CpGIs.GRanges))])

	overlapped.pos<-tapply(as.character(start(TSS))[subjectHits(overla)],queryHits(overla),paste,collapse=', ')

	DM.CpGIs.stat<-cbind(DM.CpGIs.stat,'overlapped.pos'=overlapped.pos[as.character(1:length(DM.CpGIs.GRanges))])

	ovrl.dir<-tapply(as.character(strand(TSS))[subjectHits(overla)],queryHits(overla),paste,collapse=', ')

	DM.CpGIs.stat<-cbind(DM.CpGIs.stat,'ovrl.dir'=ovrl.dir[as.character(1:length(DM.CpGIs.GRanges))])

	message('mapped')

	DM.Genes
}


