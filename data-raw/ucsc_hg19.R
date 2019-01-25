library(org.Hs.eg.db)

suppressMessages(
  geneSymbols.by.ENTEZId <- select(
  org.Hs.eg.db,
  keys=keys(org.Hs.eg.db,keytype = 'ENTREZID'),
  columns=c('SYMBOL'),
  keytype='ENTREZID'
  )
)

rownames(geneSymbols.by.ENTEZId)=geneSymbols.by.ENTEZId[,1]



