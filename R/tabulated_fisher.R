#Differential.Coverage library
#A. Favorov, E. Fertig, S. Wheelan 2014
#tabulated Fisher

#parameters: folder (default .)
#Y
#N
#theese two are numbers of positive and negative samples.
#we form the file as tabulated.fisher.Y.N.dat
#the dataframe tabulated.fisher
#each line corresponds for 4-pole table:
# MY    MN 
#MY-MY N-MN
#MY is methylated-positive  and MN is methylated-negative
#rownumber = N * MY + MN + 1 (standard 2-D transform)

prepare.tabulated.fisher<-function(Y,N,folder='.',load=TRUE,save=TRUE)
{
	while ('/'==substr(folder,nchar(folder),nchar(folder)))
		folder<-substr(folder,1,nchar(folder)-1)
	slash<-ifelse(''==folder,'','/')
	Rda.name<-paste0(folder,slash,'tabulated.fisher.',Y,'.',N,'.Rda')
	# we can the whole thing to noodles.M.Rda
	if(load && file.exists(Rda.name))
	{
		loaded<-load(Rda.name)
		if ('tabulated.fisher' %in% loaded) 
			if (class(tabulated.fisher)=='data.frame')
			{
				message('fisher tabulation loaded')
				return (tabulated.fisher)
			}
	}

	message('calcualting fisher tabulation')
	different.tests.number<-(Y+1)*(N+1)

	tabulated.fisher<-data.frame('fisher.p.values'=numeric(different.tests.number),'meth.in.normals.ratio'=numeric(different.tests.number),'meth.in.tumors.ratio'=numeric(different.tests.number),
		'OR'=numeric(different.tests.number),'CI_95_L'=numeric(different.tests.number),'CI_95_H'=numeric(different.tests.number))

	for (MY in 0:Y)
		for (MN in 0:N)
		{
			cotable<-matrix(c(MY,Y-MY,MN,N-MN),ncol=2)
			fisherres<-fisher.test(cotable)
			tabulated.fisher[N*MY+MN+1,]<-c(fisherres$p.value,cotable[2,2]/cotable[1,2],cotable[2,1]/cotable[1,1],fisherres$estimate,fisherres$conf.int[1],fisherres$conf.int[2])
		}	
	#print(dim(tabulated.fisher))
	message('saving fisher tabulation')
	if(save)
		save(file=Rda.name,list=c('tabulated.fisher','Y','N'))
		
	return (tabulated.fisher) 	
}

