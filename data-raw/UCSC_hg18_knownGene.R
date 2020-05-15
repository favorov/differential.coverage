library(TxDb.Hsapiens.UCSC.hg18.knownGene)
library(rtracklayer)
library(dplyr)
library(usethis)
library(org.Hs.eg.db)

suppressMessages(
	geneSymbols.by.ENTREZId <- AnnotationDbi::select(
	org.Hs.eg.db,
	keys=keys(org.Hs.eg.db,keytype = 'ENTREZID'),
	columns=c('SYMBOL'),
	keytype='ENTREZID'
	)
)

#could be: 
#		genelist<-unlist(genes(
#				get(.USCS.knownGenes.by.genome.annotation.id(genome.annotation.id)),
#				single.strand.genes.only=FALSE


genelist<-genes(
				TxDb.Hsapiens.UCSC.hg18.knownGene,
				single.strand.genes.only=TRUE
			)

#remove genes on strange chromosomes 
genelist<-genelist[nchar(as.character(seqnames(genelist)))<6]
#remove strange chromosomes 
seqs<-seqinfo(genelist)
seqnms<-seqnames(seqs) #all names
seqnms<-seqnms[nchar(seqnms)<6] #only real sequence names
seqs<-seqs[seqnms] #seqinfo is subesettable only by names, we make useful chromosome list
#we neeeded to restore gene_id field for each gene
genelist$gene_id=names(genelist)
#we return GRanges
ucsc_hg18_genes<-GRanges(
	ranges = ranges(genelist),
	seqnames = as.character(seqnames(genelist)),strand=strand(genelist),
	seqinfo=seqs,
	gene_id=genelist$gene_id,
	gene_name=geneSymbols.by.ENTREZId[genelist$gene_id,2]
)

usethis::use_data(ucsc_hg18_genes,overwrite=TRUE)
