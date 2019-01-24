#differential.coverage library
#A. Favorov, E. Fertig, D.Gaykalova, J. Califano, S. Wheelan 2014-2019
#chromosom features utilities

#'nucl.chromosomes.hg19
#'
#'Returns the choromosome lengths for hg19
#'
#'Returns the choromosome lengths for hg19 as a \code{seqinfo} object.
#'@export
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

#'nucl.chromosomes.hg18
#'
#'Returns the choromosome lengths for hg18
#'
#'Returns the choromosome lengths for hg18 as a \code{seqinfo} object.
#'@export
#'@param chrM boolean: if FALSE (default), returns info only for chr 1..22, X and Y. If TRUE, 1..22, X, Y and M
#'@return \code{seqinfo} object with the chromosome desriptions
#'@seealso seqinfo
nucl.chromosomes.hg18<-function(chrM=FALSE)
{
	chr.all<-seqinfo(TxDb.Hsapiens.UCSC.hg18.knownGene)
	if (chrM)
		chr.all[names(chr.all)[1:25]]
	else
		chr.all[names(chr.all)[1:24]]
}

#'nucl.chromosomes.hg38
#'
#'Returns the choromosome lengths for hg38
#'
#'Returns the choromosome lengths for hg18 as a \code{seqinfo} object.
#'@export
#'@param chrM boolean: if FALSE (default), returns info only for chr 1..22, X and Y. If TRUE, 1..22, X, Y and M
#'@return \code{seqinfo} object with the chromosome desriptions
#'@seealso seqinfo
nucl.chromosomes.hg38<-function(chrM=FALSE)
{
	chr.all<-seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene)
	if (chrM)
		chr.all[names(chr.all)[1:25]]
	else
		chr.all[names(chr.all)[1:24]]
}


#'get.cytoband.ranges for hg19, no hg18 or hg 38 version yet
#'
#'Returns \code{GRanges} with cytobands.
#'
#'Returns \code{GRanges} with cytobands. Actually, takes \code{\link[SNPchip]{getCytoband}}  from \pkg{SNPchip} and convert it into \code{GRanges}. \code{getCytoband} returns bands in UCSC notation (0-based start, 1-based end), get.cytoband.ranges return it in 1-based.
#'
#'@export
#'@seealso SNPchip, getCytoband 
get.cytoband.ranges<-function()
{
	cbands<-getCytoband()
	karyotype<-GRanges(
		seqinfo=nucl.chromosomes.hg19(),
		seqnames=cbands$chrom,
		range=IRanges(start=cbands$start+1,end=cbands$end),
		name=cbands$name,
		gieStain=cbands$gieStain)
	karyotype[order(karyotype)]
}
