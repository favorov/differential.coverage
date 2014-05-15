#Differential.Coverage library
#A. Favorov, S. Wheelan 2014
#annotation utilities

nucl.chromosomes.hg19<-function()
{
	chr.all<-seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)
	chr.all[names(chr.all[1:24])]
	#only chrNN, chrX and chrY.
	#chrM is 25
}
