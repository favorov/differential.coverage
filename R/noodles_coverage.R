#Differential.Coverage library
#A. Favorov, E. Fertig, S. Wheelan 2014

#'CountCoverageOfNoodles
#'
#'It is a central function of all the Differential.Coverage package. It gets a set of intervals (noodles) as a \code{GRanges} and a list of names of bed files.
#'Each bedfile represents a sample. For each noodle and each sample, the lenght of covered part of the noodle is calculated.
#'
#'@param noodles \code{GRanges} with the intervals
#'@param bedfilnames list of names of bedfiles, one per sample, with some (e.g. methylation) coverage information
#'@param bed.ids optional list of names for the samples, they will be used as column names in the result. Default = \code{bedfilnames}
#'@return \code{data.frame}, each row correspond to a noodle; columns are samples
#'@seealso \code{Differential.Coverage}
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
