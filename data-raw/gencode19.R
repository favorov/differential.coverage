library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(rtracklayer)

if(! 'raw_differtial_coverage_gencode_19_gff' %in% ls()){
  raw_differtial_coverage_gencode_19_gff=
      readGFF("http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz",
              version=3)
}
