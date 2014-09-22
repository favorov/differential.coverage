#differential.coverage library
#A. Favorov, E. Fertig, S. Wheelan 2014
#prepare noodles

#returns noodles that uniformly cover
#the set of chromosomes that are described by 
#sequinfo
#seqinfo is seqinfo or a vector of ints with the names as chromosome.names
prepare.uniform.noodles<-function(seq.info,noodle.length=1000)
{
	
	if (class(seq.info) == 'Seqinfo')
		seq.info<-seqlengths(seq.info)
	#now, seqinfo is a vector of ints with the names
	noodles.seqnames=character(0)
	noodles.start=integer(0)
	noodles.end=integer(0)
	message('generating noodles')
	for (chr in names(seq.info))
	{
		this.noodles.start<-seq(1,seq.info[chr],by=noodle.length)
		noodles.start<-c(noodles.start,this.noodles.start)
		
		noodles.seqnames<-c(noodles.seqnames,rep(chr,length(this.noodles.start)))
		
		this.noodles.end<-this.noodles.start+noodle.length-1
		this.noodles.end[length(this.noodles.end)]<-
			min(this.noodles.end[length(this.noodles.end)],seq.info[chr])
		#for the last noodle not to break end of chr
		noodles.end<-c(noodles.end,this.noodles.end)
	}
	message('combining noodles')
	GRanges(
		seqlengths=seq.info,
		seqnames=noodles.seqnames,
		ranges=IRanges
		(
			start=noodles.start,
			end=noodles.end
		)
	)
}
