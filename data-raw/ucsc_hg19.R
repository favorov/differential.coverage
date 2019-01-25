library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(usethis)

suppressMessages(
  geneSymbols.by.ENTEZId <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys=keys(org.Hs.eg.db,keytype = 'ENTREZID'),
  columns=c('SYMBOL'),
  keytype='ENTREZID'
  )
)

rownames(geneSymbols.by.ENTEZId)=geneSymbols.by.ENTEZId[,1]

genes<-unlist(genes(TxDb.Hsapiens.UCSC.hg19.knownGene,single.strand.genes.only=FALSE))
#remove strange threads
genes<-genes[nchar(as.character(seqnames(genes)))<6]






