library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(rtracklayer)
library(dplyr)
library(usethis)

#gencode rel 29, hg38
#https://www.gencodegenes.org/human/release_29.html

if(! 'raw_differtial_coverage_gencode_hs29_gff' %in% ls()){
  raw_differtial_coverage_gencode_hs29_gff=
      readGFF("http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gff3.gz",
              version=3)
}

# we use UCSC seqinfo from TxDb.Hsapiens.UCSC.hg38.knownGene
seqs<-seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene)
seqnms<-seqnames(seqs) #all names
seqnms<-seqnms[nchar(seqnms)<6] #only real sequence names
seqs<-seqs[seqnms] #seqinfo is subesettable only by names, we make useful chromosome list
gencode_hs29_genes<-
  GRanges(as.data.frame(raw_differtial_coverage_gencode_hs29_gff) %>% 
      filter(type=="gene") %>% 
      select(seqid,start,end,strand,gene_name,ensembl=gene_id),
  seqinfo = seqs)
#gene_id here is ensembl name, but withe additional .n (transcript#), we remove it
gencode_hs29_genes$ensembl<-sapply(strsplit(gencode_hs29_genes$ensembl,'.',fixed = TRUE),"[[",1)

#we are in the package root folder
usethis::use_data(gencode_hs29_genes,overwrite = TRUE)
