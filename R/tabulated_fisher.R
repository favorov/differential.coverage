#Differential.Coverage library
#A. Favorov, E. Fertig, S. Wheelan 2014
#tabulated Fisher

#'prepare.tabulated.fisher
#'
#'forms the dataframe (\code{tabulated.fisher}) with the tabulted results for all possible Fisshe exact tests with Y and N numbers of cases and controls
#'
#'each row gives statistics (p-value, etc) for a 4-pole table:
#'\tabular{cc}{
#'MY \tab MN \cr 
#'Y-MY \tab N-MN }
#'MY is methylated-positive and MN is methylated-negative (normal)
#'we want to represent all the Fishers as a (Y+1) rows * (N+1) cols
#'matrix.
#'MY=0..Y
#'MN=0..N
#'then, we map it to 1..(MY+1)(MN+1), and it is the row index
#'@param Y cases# , e.g. number of cancer samples 
#'@param N controls#
#'@param folder (default .) 
#'@param load (default TRUE) logical. If TRUE, we try to load file tabulated.fisher.Y.N.dat
#'@param save (default TRUE) logical. If TRUE, we save the \code{tabulated.fisher} to file tabulated.fisher.Y.N.dat
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

	#again: Y=#of positive (tumors)
	#N=#of negative 
	#MY MN
	#Y-MY N-MN
	for (MY in 0:Y)
		for (MN in 0:N)
		{
			cotable<-matrix(c(MY,MN,Y-MY,N-MN),ncol=2,byrow=TRUE)
			fisherres<-fisher.test(cotable)
			tabulated.fisher[tab.fisher.row.no(Y,N,MY,MN),]<-c(fisherres$p.value,cotable[1,2]/N,cotable[1,1]/Y,fisherres$estimate,fisherres$conf.int[1],fisherres$conf.int[2])
		}	
	#print(dim(tabulated.fisher))
	message('saving fisher tabulation')
	if(save)
		save(file=Rda.name,list=c('tabulated.fisher','Y','N'))
	return (tabulated.fisher) 	
}


#'tab.fisher.row.no
#'
#'calculates the index of the row (each row is a Fisher table)
#'that represents the Fisher's table among the matrix that tabulate the Fisher table results
#'
#'@param Y sum of column 1 (cases#) 
#'@param N sum of column 2 (controls#)
#'@param MY methylated cases #
#'@param MN methylated controls #

tab.fisher.row.no<-function(Y,N,MY,MN)
{
	(N+1)*MY+MN+1	
}

