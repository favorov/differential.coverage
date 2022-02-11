library(GenomeInfoDb)
library(rtracklayer)
library(dplyr)
library(usethis)

#gencode rel 23, current for mm10 (GRCm38)
#https://www.gencodegenes.org/mouse/release_M23.html

if(! 'raw_differtial_coverage_gencode_mm23_gff' %in% ls()){
  raw_differtial_coverage_gencode_mm32_gff=
      readGFF("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gff3.gz",
              version=3)
}

# we use UCSC seqinfo from TxDb.Hsapiens.UCSC.hg38.knownGene
seqs<-Seqinfo(genome="mm10")
seqnms<-seqnames(seqs) #all names
seqnms<-seqnms[nchar(seqnms)<6] #only real sequence names
seqs<-seqs[seqnms] #seqinfo is subesettable only by names, we make useful chromosome list
gencode_mm23_genes<-
  GRanges(as.data.frame(raw_differtial_coverage_gencode_mm23_gff) %>% 
      filter(type=="gene") %>% 
      select(seqid,start,end,strand,gene_name,ensembl=gene_id),
  seqinfo = seqs)
#gene_id here is ensembl name, but withe additional .n (transcript#), we remove it
gencode_mm23_genes$ensembl<-sapply(strsplit(gencode_mm23_genes$ensembl,'.',fixed = TRUE),"[[",1)

#we are in the package root folder
usethis::use_data(gencode_mm23_genes,overwrite = TRUE)
