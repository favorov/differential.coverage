#Differential.Coverage library
#A. Favorov, S. Wheelan 2014

CountCoverageOfNoodles<-function(noodles,bedfilnames,bed.ids=bedfilnames){
	if (length(bedfilnames)!=length(bed.ids))
	{
		stop("The lists of bed file names and bed ids has different lenghts. So what?")
	}
	else
		names(bedfilnames)=bed.ids
	for (bed.id in bed.ids)
	{
		print(bed.id)
	}
}
