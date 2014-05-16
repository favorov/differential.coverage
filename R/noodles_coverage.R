#Differential.Coverage library
#A. Favorov, S. Wheelan 2014

CountCoverageOfNoodles<-function(noodles,bedfilnames,bed.ids=bedfilnames){
	if (class(noodles)!='GRanges')
	{
		stop("The noodles for coverage is not GRanges. So what?")
	}
	if (length(bedfilnames)!=length(bed.ids))
	{
		stop("The lists of bed file names and bed ids has different lenghts. So what?")
	}
	else
		names(bedfilnames)<-bed.ids
	message('coverage')
	noodles.coverage<-data.frame(matrix(0,ncol=length(bed.ids),nrow=length(noodles)))
	colnames(noodles.coverage)<-bed.ids
	for (bed.id in bed.ids)
	{
		message(bed.id)
		beads<-import(bedfilnames[bed.id])
		overrle<-findOverlaps(noodles,beads)
		covered<-tapply(width(beads[subjectHits(overrle)]),queryHits(overrle),sum)
		noodles.coverage[as.integer(names(covered)),bed.id]<-covered
	}
	noodles.coverage
}
#it counts and return noodles.coverage data frame
