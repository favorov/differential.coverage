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
#rownumber = (N-1) MY + MN (standard 2-D transform)

prepare.tabulated.fisher(Y,N)
{
	Rda.name<-paste0('tabulated.fisher.',Y,'.',N,'.Rda')
}
