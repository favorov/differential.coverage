#differential.coverage library
#A. Favorov, E. Fertig, D.Gaykalova, J. Califano, S. Wheelan 2014-2020
#annotation utilities

#' @import rtracklayer
#' @import stringr
#import description end 
0

#' get.Known.Gene.List
#' 
#' refer to a known gene list for the annotation functions
#' corrently, we have: gencode hs 19 (current for hg 19), gencode hs 26, gencode hs 29, gencode hs 32 (current for hg38), gencode hs 34; 
#' ucsc hg 18; uscs hg 19; ucsc hg 38 (all in single-strand option)
#' 
#' @export
#' 
#' @param genome.annotation.id says what annotation is used to prepare the output. \code{gencode19} (default), \code{gencode29} and \code{gencode26} load genes for gencode stable annotations (19 and 29) and for 26 (gencode 19 is for human genome 19, gencode 26 and 29 are for 38 version). \code{hg18}, \code{hg19}, and \code{hg38} load \code{TxDb.Hsapiens.UCSC.hg18.knownGene}, \code{TxDb.Hsapiens.UCSC.hg19.knownGene} and \code{TxDb.Hsapiens.UCSC.hg38.knownGene}, correspondingly. 
##' @param single.strand.genes.only UCSC annotations contain ~500 pair of same-named genes that exist on both strands.The parameter says whether to exclude them from the gene list to be returned. The default is \code{FALSE} that allows these genes to be included.
#' @return \code{GRanges} object that contains the gene annotation. The gene_name metadata field is the gene symbol according to the requested annotation. ucsc* uses org.Hs.eg.db names, the gencode provides its own gene names 
get.Known.Gene.List<-function(genome.annotation.id='gencode34')
#single.strand.genes.only=FALSE)
{
	human<-FALSE
	mouse<-FALSE
	ucsc<-FALSE
	gencode<-FALSE
	number<-0
	#defaults
	genome_id<-str_to_lower(genome.annotation.id)
	if ( str_detect(genome_id,'hs') || str_detect(genome_id,'hg') || str_detect(genome_id,'human')) human<-TRUE
	if ( str_detect(genome_id,'mm') || str_detect(genome_id,'mm') || str_detect(genome_id,'mouse')) mouse<-TRUE
	if ( str_detect(genome_id,'ucsc') ) ucsc<-TRUE
	if ( str_detect(genome_id,'gencode') ) gencode<-TRUE 
	digi<-str_extract_all(genome_id,"[:digit:]+")[[1]] #only one string is parsed
	if(length(digi) > 1 ) stop("More that one number in genome id. So what?")
	digi<-digi[1] #the only number
	digi<-str_c("n",digi) #n28" instead of "28"
	genome<-switch(digi,
		n9={if (human || gencode) NULL else ucsc_mm9_genes},
		n10={if (human || gencode) NULL else ucsc_mm10_genes},
		n23={if (human || ucsc) NULL else gencode_mm23_genes},
		n24={if (human || ucsc) NULL else gencode_mm24_genes},
		n25={if (human || ucsc) NULL else gencode_mm25_genes},
		n18={if (mouse || gencode) NULL else ucsc_hg18_genes},
		n19={if (mouse || ! (gencode || ucsc)) stop (str_c("We are not sure whether you nean hg19 ucsc annotation or hs19 gencode : \"",genome_id,"\"."))	else if (ucsc) ucsc_hg19_genes else gencode_hs19_genes},
		n38={if (mouse || gencode) NULL else ucsc_hg38_genes},
		n26={if (mouse || ucsc) NULL else gencode_hs26_genes},
		n29={if (mouse || ucsc) NULL else gencode_hs29_genes},
		n32={if (mouse || ucsc) NULL else gencode_hs32_genes},
		n34={if (mouse || ucsc) NULL else gencode_hs34_genes}
	)
	if (is.null(genome )) stop (str_c("Unknown genome : \"",genome_id,"\"."));
	genome
}

#' inflate.noodles
#' 
#' Expands each element of GRanges by flanks value both sides. Check chomosome boundaries and avoid breaking them (thus differs from \link{resize}). Needs seqlengths for chromosemes to be defined in GRanges or by seqlengths parameter. If both are given, check whether they do not conradict.  
#' 
#' @export
#' 
#' @param noodles the \code{GRanges} list of intervals to inflate
#' @param flanks lenght to inflate the noddles by before the search; if >0, the seqlengths information is to be set as \code{seqlengths} or in \code{seqlengths(noodles)}
#' @param seqlengths is to provide the chromosome lentth information without including it in noodles
#' 
inflate.noodles<-function
(
	noodles, # GRanges with the noodles
	flanks=0, # how far to inflate 
	seqlengths=NA #chromosome length, to be subsettable by chomosome name as charater if not NA
)
{
	inflated.noodles<-noodles

	if (sum(is.na(seqlengths))) 
	#we did not provide the lengths explicitely
	#sum is to supress warnong is seqlengths is provided as vector
	#for scalar NA sum is the same as its agrument
	{
		seqlengths<-seqlengths(noodles)
	}
	if (flanks>0) #valudate seqlengths
	{
		if(sum(is.na(seqlengths[as.character(seqnames(noodles))]))>0)
			stop('Inflating noodles, seqlengths cannot be undefined here')
	}	
	#inflate noodles
	start(inflated.noodles)<-pmax(1,start(noodles)-flanks)
	end(inflated.noodles)<-pmin(end(noodles)+flanks,as.integer(seqlengths[as.character(seqnames(noodles))]))
	inflated.noodles
}



#' genes.with.TSS.covered
#' 
#' Generates list of genes that start inside a given set of intervals
#' 
#' After the noodles (the set of intervals to search TSS in) are inflated by flanks, we look for all the TSS that start inside the (inflated) intervals. Genes are provided by parameter or they are returned by \code{\link{get.Known.Gene.List}} for \code{genome_id}. If the noodles has p.value and/or fdr metadata, we ascribe the data of the interval to the retrieved gene. If there a gene refers to a set of noodles, it has min ascribed. The ishyper data is also transferred to gene, if it is not contradictory.
#' 
#' @export
#' 
#' @param noodles the \code{GRanges} list of intervals to look TSS in
#' @param flanks lenght to inflate the noddles by before the search; if >0, the seqlengths information is to be set in \code{noodles}
#' @param genes it is GRandes oject with the annotation genes. The gene_name matadata field is almost required, it passes the gene name to form output. 
#' The default is NA that means that the function will call \code{\link{get.Known.Gene.List}} 
#' @param genome.id is genome.annotation.id to call \code{\link{get.Known.Gene.List}} 
#' @return \code{GRanges} object that is the list of the genes we look for - the object is not co-indexed with \code{noodles} parameter
#' @author Alexander Favorov, \email{favorov@@sensi.org}
#' @seealso \code{\link{differential.coverage}}
#' @keywords annotation
genes.with.TSS.covered<-function(
	noodles, # GRanges with the noodles, if it has p.value, fdr and ishyper values, they will be mapped to genes
	flanks=0, #how far to shrink
	genes=NA,
	genome.id='gencode19'
)
{
	#prepare genelist
	if(!is.na(genes) && class(genes)=='GRanges') {genelist<-genes} 
	else{genelist<-get.Known.Gene.List(genome.id)}
		
	#inflate noodles
	inflated.noodles<-inflate.noodles(noodles,flanks,seqlengths(genelist))

	#initialise the list to subset later
	#the list is a Granges
	TSS<-genelist

	tss.start<-ifelse(strand(TSS)=='+',start(TSS),end(TSS))
	start(TSS)<-tss.start
	end(TSS)<-tss.start
	#TSS prepared
	message("tss prepared")

	#make overlap
	overlaps<-findOverlaps(inflated.noodles,TSS)
	message('overlapped')

	noodle.TSS.Gene.Indices<-unique(subjectHits(overlaps))
	#indices of all the genes tha were hit by any noodle
	#maybe, some of them are hit with more than one

	noodle.TSS.Genes<-genelist[noodle.TSS.Gene.Indices]
	

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

#' genes.intersected
#' 
#' Generates list of genes that intersect a given set of intervals
#' 
#' After the noodles (the set of intervals to search intersection with) are inflated by flanks, we look for all the genes that intersect the (inflated) intervals. Genes are provided by parameter or they are returned by \code{\link{get.Known.Gene.List}} for \code{genome_id}. If the noodles has p.value and/or fdr metadata, we ascribe the data of the interval to the retrieved gene. If there a gene refers to a set of noodles, it has min ascribed. The ishyper data is also transferred to gene, if it is not contradictory.
#' @export
#' @inheritParams genes.with.TSS.covered
#' @return \code{GRanges} object that is the list of the genes we look for - the object is not co-indexed with \code{noodles} parameter
#' @inherit genes.with.TSS.covered author seealso 
#' @keywords annotation
genes.intersected<-function(
	noodles, # GRanges with the noodles, if it has p.value, fdr and ishyper values, they will be mapped to genes
	flanks=0, #how far to shrink
	genes=NA,
	genome.id='gencode19'
)
{
	#prepare genelist
	if(!is.na(genes) && class(genes)=='GRanges') {genelist<-genes} 
	else{genelist<-get.Known.Gene.List(genome.id)}
	
	#inflate noodles
	inflated.noodles<-inflate.noodles(noodles,flanks,seqlengths(genelist))

	#make overlap
	overlaps<-findOverlaps(inflated.noodles,genelist)
	message('overlapped')

	noodle.Gene.Indices<-unique(subjectHits(overlaps))
	#indices of all the genes tha were hit by any noodle
	#maybe, some of them are hit with more than one

	noodle.Genes<-genelist[noodle.Gene.Indices]
	
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


#' genes.with.TSS.covered.by.interval
#' 
#' Generates a list of genes (possibly, empty) that start inside each interval
#' 
#' After the noodles (the set of intervals to search TSS in) are inflated by flanks, we look for all the TSS that start inside the (inflated)  intervals. Genes are provided by parameter or they are returned by \code{\link{get.Known.Gene.List}} for \code{genome_id}. If a noodle (interval) overlaps more that one TSS, we form a text list of the genes.
#' @export
#' @inheritParams genes.with.TSS.covered
#' @return \code{GRanges} object, noodles argument with added TSS-overlapped genes for each interval  
#' @inherit genes.with.TSS.covered author seealso 
#' @keywords annotation
genes.with.TSS.covered.by.interval<-function(
	noodles, # GRanges with the noodles, if it has p.value ans ishyper values, thay will be mapped to genes
	flanks=0, #how far to shrink 
	genes=NA,
	genome.id='gencode19'
)
{
	#prepare genelist
	if(!is.na(genes) && class(genes)=='GRanges') {genelist<-genes} 
	else{genelist<-get.Known.Gene.List(genome.id)}
		
	#inflate noodles
	inflated.noodles<-inflate.noodles(noodles,flanks,seqlengths(genelist))
	
	#initialise the list to subset later
	#the list is a Granges
	TSS<-genelist
	
	tss.start<-ifelse(strand(TSS)=='+',start(TSS),end(TSS))

	start(TSS)<-tss.start
	end(TSS)<-tss.start
	#IAMHERE
	#make overlap
	overlapa<-findOverlaps(inflated.noodles,TSS)
	message('overlapped')
	
	decorated.noodles<-noodles

	#overlapped.TSS

	overlapped.TSS<-tapply(TSS$gene_name[subjectHits(overlapa)],queryHits(overlapa),paste,collapse=', ')
	decorated.noodles$overlapped.TSS<-overlapped.TSS[as.character(1:length(decorated.noodles))]
	#we need this addressing scheme ([as.character(1:length(decorated.noodles))]) because names(overlapped.TSS)
	#are the indices of noodles that have overlapped TSS. Those that have not are not represented in overlapped.TSS
	#and after our indexing they are NA, and it is exacltly what we want them to be

	overlapped.pos<-tapply(as.character(start(TSS))[subjectHits(overlapa)],queryHits(overlapa),paste,collapse=', ')
	decorated.noodles$overlapped.pos<-overlapped.pos[as.character(1:length(decorated.noodles))]

	ovrl.dir<-tapply(as.character(strand(TSS))[subjectHits(overlapa)],queryHits(overlapa),paste,collapse=', ')
	decorated.noodles$ovrl.dir<-ovrl.dir[as.character(1:length(decorated.noodles))]
	
	if ("ensembl" %in% colnames(mcols(TSS))) {
	  overlapped.TSS.ensembl<-tapply(TSS$ensembl[subjectHits(overlapa)],queryHits(overlapa),paste,collapse=', ')
	  decorated.noodles$overlapped.TSS.ensembl<-overlapped.TSS.ensembl[as.character(1:length(decorated.noodles))]
	}
	
	if ("gene_id" %in% colnames(mcols(TSS))) {
	  overlapped.TSS.gene_id<-tapply(TSS$gene_id[subjectHits(overlapa)],queryHits(overlapa),paste,collapse=', ')
	  decorated.noodles$overlapped.TSS.gene_id<-overlapped.TSS.gene_id[as.character(1:length(decorated.noodles))]
	}
	
	message('mapped')

	decorated.noodles

}


#' genes.intersected.by.interval
#' 
#' Generates a list of genes (possibly, empty) that intersects with each interval
#' 
#' After the noodles (the set of intervals to search intersections with) are inflated by flanks, we look for all the genes that overlap with the (inflated) intervals. Genes are provided by parameter or they are returned by \code{\link{get.Known.Gene.List}} for \code{genome_id}. If a noodle (interval) overlaps more that one gene, we form a text list of the genes.
#' @export
#' @inheritParams genes.with.TSS.covered
#' @return \code{GRanges} object, noodles argument with added gene-overlapped genes for each interval  
#' @inherit genes.with.TSS.covered author seealso 
#' @keywords annotation
genes.intersected.by.interval<-function(
	noodles, # GRanges with the noodles, if it has p.value ans ishyper values, thay will be mapped to genes
	flanks=0, #how far to shrink 
	genes=NA,
	genome.id='gencode19'
)
{
	#prepare genelist
	if(!is.na(genes) && class(genes)=='GRanges') {genelist<-genes} 
	else{genelist<-get.Known.Gene.List(genome.id)}
		
	#inflate noodles
	inflated.noodles<-inflate.noodles(noodles,flanks,seqlengths(genelist))
	
	
	overlapa<-findOverlaps(inflated.noodles,genelist)
	message('overlapped')
	
	decorated.noodles<-noodles

	#overlapped.genes

	overlapped.genes<-tapply(genelist$gene_name[subjectHits(overlapa)],queryHits(overlapa),paste,collapse=', ')
	decorated.noodles$overlapped.genes=overlapped.genes[as.character(1:length(decorated.noodles))]
	#we need this addressing scheme ([as.character(1:length(decorated.noodles))]) because names(overlapped.TSS)
	#are the indices of noodles that have overlapped TSS. Those that have not are not represented in overlapped.TSS
	#and after our indexing they are NA, and it is exacltly what we want them to be

	overlapped.pos<-tapply(as.character(start(genelist))[subjectHits(overlapa)],queryHits(overlapa),paste,collapse=', ')
	decorated.noodles$overlapped.pos<-overlapped.pos[as.character(1:length(decorated.noodles))]

	ovrl.dir<-tapply(as.character(strand(genelist))[subjectHits(overlapa)],queryHits(overlapa),paste,collapse=', ')
	decorated.noodles$ovrl.dir<-ovrl.dir[as.character(1:length(decorated.noodles))]

	if ("ensembl" %in% colnames(mcols(genelist))) {
	  overlapped.genelist.ensembl<-tapply(genelist$ensembl[subjectHits(overlapa)],queryHits(overlapa),paste,collapse=', ')
	  decorated.noodles$overlapped.genelist.ensembl<-overlapped.genelist.ensembl[as.character(1:length(decorated.noodles))]
	}
	
	if ("gene_id" %in% colnames(mcols(genelist))) {
	  overlapped.genelist.gene_id<-tapply(genelist$gene_id[subjectHits(overlapa)],queryHits(overlapa),paste,collapse=', ')
	  decorated.noodles$overlapped.genelist.gene_id<-overlapped.genelist.gene_id[as.character(1:length(decorated.noodles))]
	}
	message('mapped')

	decorated.noodles

}


#' closest.gene.by.interval
#' 
#' Find the closest gene for each of the given set of intervals.  Genes are provided by parameter or they are returned by \code{\link{get.Known.Gene.List}} for \code{genome_id}. 
#' 
#' @export
#' @inheritParams genes.with.TSS.covered
#' @return \code{GRanges} object, noodles argument with added closest gene info for each interval  
#' @inherit genes.with.TSS.covered author seealso 
#' @keywords annotation
closest.gene.by.interval<-function(
	noodles, # GRanges with the noodles, if it has p.value ans ishyper values, thay will be mapped to genes
	genes=NA,
	genome.id='gencode19'
)
{
	#prepare genelist
	if(!is.na(genes) && class(genes)=='GRanges') {genelist<-genes} 
	else{genelist<-get.Known.Gene.List(genome.id)}
		
	message('closest')

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

	noodles.decoration$closest.gene[is.a.near.gene]<-genelist$gene_name[near.gene]
	noodles.decoration$start[is.a.near.gene]<-start(genelist)[near.gene]
	noodles.decoration$end[is.a.near.gene]<-end(genelist)[near.gene]
	noodles.decoration$strand[is.a.near.gene]<-
		as.character(strand(genelist)[near.gene])
	noodles.decoration$dist[is.a.near.gene]<-dist.gene


	if ("ensembl" %in% colnames(mcols(genelist))) {
	  noodles.decoration$ensembl<-character(length(noodles))
		noodles.decoration$ensembl[is.a.near.gene]<-genelist$ensembl[near.gene]
	}
	
	if ("gene_id" %in% colnames(mcols(genelist))) {
	  noodles.decoration$gene_id<-character(length(noodles))
		noodles.decoration$gene_id[is.a.near.gene]<-genelist$gene_id[near.gene]
	}
	
	noodles.decoration[!is.a.near.gene,]=NA

	decorated.noodles<-noodles

	mcols(decorated.noodles)<-noodles.decoration

	message('mapped')
	
	decorated.noodles
}


#' closest.gene.start.by.interval
#' 
#' Find the closest gene start for each of the given set of intervals.  Genes are provided by parameter or they are returned by \code{\link{get.Known.Gene.List}} for \code{genome_id}. 
#' 
#' @export
#'
#' @inheritParams genes.with.TSS.covered
#' @return \code{GRanges} object, noodles argument with added closest gene info for each interval  
#' @inherit genes.with.TSS.covered author seealso 
#' @keywords annotation
closest.gene.start.by.interval<-function(
	noodles, # GRanges with the noodles, if it has p.value ans ishyper values, thay will be mapped to genes
	genes=NA,
	genome.id='gencode19'
)
{
	#prepare genelist (we call it TSS here)
	if(!is.na(genes) && class(genes)=='GRanges') {TSS<-genes} 
	else{TSS<-get.Known.Gene.List(genome.id)}
	
	message('closest')

	#prepare gene TSS; we refere the TxDb object by name

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

	noodles.decoration$closest.TSS[is.a.near.TSS]<-TSS$gene_name[near.TSS]
	noodles.decoration$pos[is.a.near.TSS]<-start(TSS)[near.TSS]
	noodles.decoration$dir[is.a.near.TSS]<-as.character(strand(TSS)[near.TSS])
	noodles.decoration$dist[is.a.near.TSS]<-dist.TSS

	qqTSS<<-TSS
	
	if ("ensembl" %in% colnames(mcols(TSS))) {
	  noodles.decoration$ensembl<-character(length(noodles))
		noodles.decoration$ensembl[is.a.near.TSS]<-TSS$ensembl[near.TSS]
	}
	
	if ("gene_id" %in% colnames(mcols(TSS))) {
	  noodles.decoration$gene_id<-character(length(noodles))
		noodles.decoration$gene_id[is.a.near.TSS]<-TSS$gene_id[near.TSS]
	}
	
	noodles.decoration[!is.a.near.TSS,]=NA

	decorated.noodles<-noodles

	mcols(decorated.noodles)<-noodles.decoration

	message('mapped')
	
	decorated.noodles
}


#' preceded.gene.by.interval
#' 
#' For each interval, find the closest gene that the interval precede.  
#' An overlapped gene is neglected.
#' Genes are provided by parameter or they are returned by \code{\link{get.Known.Gene.List}} for \code{genome_id}. 
#' 
#' @export
#'
#' @inheritParams genes.with.TSS.covered
#' @return \code{GRanges} object, noodles argument with added closest preceded (non overlapped) gene info for each interval  
#' @inherit genes.with.TSS.covered author seealso 
#' @keywords annotation
preceded.gene.by.interval<-function(
	noodles, # GRanges with the noodles, if it has p.value ans ishyper values, thay will be mapped to genes
	genes=NA,
	genome.id='gencode19'
)
{
	#prepare genelist
	if(!is.na(genes) && class(genes)=='GRanges') {genelist<-genes} 
	else{genelist<-get.Known.Gene.List(genome.id)}
		
	message('closest preceded')

	near.gene<-precede(noodles,genelist)

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

	noodles.decoration$closest.gene[is.a.near.gene]<-genelist$gene_name[near.gene]
	noodles.decoration$start[is.a.near.gene]<-start(genelist)[near.gene]
	noodles.decoration$end[is.a.near.gene]<-end(genelist)[near.gene]
	noodles.decoration$strand[is.a.near.gene]<-
		as.character(strand(genelist)[near.gene])
	noodles.decoration$dist[is.a.near.gene]<-dist.gene


	if ("ensembl" %in% colnames(mcols(genelist))) {
	  noodles.decoration$ensembl<-character(length(noodles))
		noodles.decoration$ensembl[is.a.near.gene]<-genelist$ensembl[near.gene]
	}
	
	if ("gene_id" %in% colnames(mcols(genelist))) {
	  noodles.decoration$gene_id<-character(length(noodles))
		noodles.decoration$gene_id[is.a.near.gene]<-genelist$gene_id[near.gene]
	}
	
	noodles.decoration[!is.a.near.gene,]=NA

	decorated.noodles<-noodles

	mcols(decorated.noodles)<-noodles.decoration

	message('mapped')
	
	decorated.noodles
}
