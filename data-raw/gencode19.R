library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(rtracklayer)
library(usethis)

if(! 'raw_differtial_coverage_gencode_19_gff' %in% ls()){
  raw_differtial_coverage_gencode_19_gff=
      readGFF("http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz",
              version=3)
}

if(! 'gencode19_genes' %in% ls()){
  # we use UCSC seqinfo from TxDb.Hsapiens.UCSC.hg19.knownGene
  seqs<-seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)
  seqnms<-seqnames(seqs) #all names
  seqnms<-seqnms[nchar(seqnms)<6] #only real sequence names
  seqs<-seqs[seqnms] #seqinfo is subesettable only by names, we make useful chromosome list
  gencode19_genes<-
    GRanges(as.data.frame(raw_differtial_coverage_gencode_19_gff) %>% 
        filter(type=="gene") %>% 
        select("seqid","start","end","gene_name","gene_id"),
    seqinfo = seqs)
}
#we are in the package root folder
use_data(gencode19_genes)
