#differential.coverage library
#A. Favorov, E. Fertig, D.Gaykalova, J. Califano, S. Wheelan 2014-2016

#'count.coverage.of.noodles
#'
#'It is a central function of all the differential.coverage package. It gets a set of intervals (noodles) as a \code{GRanges} and a list of names of bed files.
#'Each bedfile represents a sample. For each noodle and each sample, the the sum of the lenght of intervals (read from the sample's bed file) that intersect with the noodle is returned.
#'
#'@param noodles \code{GRanges} with the intervals
#'@param bedfilnames list of names of bedfiles, one per sample, with some (e.g. methylation) coverage information
#'@param bed.ids optional list of names for the samples, they will be used as column names in the result. Default = \code{bedfilnames}
#'@return \code{Matrix}, each row correspond to a noodle; columns are samples, sparse=TRUE
#'@seealso \code{differential.coverage}
count.coverage.of.noodles<-function(noodles,bedfilnames,bed.ids=bedfilnames){
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
	noodles.coverage<-Matrix(0,ncol=length(bed.ids),nrow=length(noodles),sparse = TRUE)
	colnames(noodles.coverage)<-bed.ids
	for (bed.id in bed.ids)
	{
		message(bed.id)
		beads<-import(bedfilnames[bed.id])
		overrle<-findOverlaps(noodles,beads)
		#end(beads)[subjectHits(overrle) is a vector of the same lenght as the hit list, and 
		#it carries the ends of all beads (actully, the bedfile intervals).
		#analogous are end(noodles)[queryHits(overrle)] (ends of all the noddles that are intersected by beads, etc
		ovelap.width<-pmin(end(beads)[subjectHits(overrle)],end(noodles)[queryHits(overrle)])-
			pmax(start(beads)[subjectHits(overrle)],start(noodles)[queryHits(overrle)])+1
		#so, overrap.width is a list of lengthes of all the noodle x bead overlaps in the same order (and length)
		#as they are listed by findOverlaps in overrle
		covered<-tapply(ovelap.width,queryHits(overrle),sum)
		noodles.coverage[as.integer(names(covered)),bed.id]<-as.integer(covered)
	}
	noodles.coverage
}

#it counts and return noodles.coverage data frame
#'indicate.any.coverage.of.noodles
#'
#'It is the alterantive central function of all the differential.coverage package. It gets a set of intervals (noodles) as a \code{GRanges} and a list of names of bed files.
#'Each bedfile represents a sample. For each noodle and each sample, 0 or 1, which indicates whether the noodle intersects with any bed interval from the sample, is returned.
#'
#'@param noodles \code{GRanges} with the intervals
#'@param bedfilnames list of names of bedfiles, one per sample, with some (e.g. methylation) coverage information
#'@param bed.ids optional list of names for the samples, they will be used as column names in the result. Default = \code{bedfilnames}
#'@return \code{Matrix}, each row correspond to a noodle; columns are samples, sparse=TRUE
#'@seealso \code{differential.coverage}
indicate.any.coverage.of.noodles<-function(noodles,bedfilnames,bed.ids=bedfilnames){
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
	noodles.coverage<-Matrix(0,ncol=length(bed.ids),nrow=length(noodles),sparse = TRUE)
	colnames(noodles.coverage)<-bed.ids
	for (bed.id in bed.ids)
	{
		message(bed.id)
		beads<-import(bedfilnames[bed.id])
		covered<-overlapsAny(noodles,beads)
		noodles.coverage[,bed.id]<-as.integer(covered)
	}
	noodles.coverage
}

CountCoverageOfNoodles<-function(noodles,bedfilnames,bed.ids=bedfilnames){
	.Deprecated('count.coverage.of.noodles')
	count.coverage.of.noodles(noodles,bedfilnames,bed.ids)
}
