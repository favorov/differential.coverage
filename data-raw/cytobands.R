library(rtracklayer)
library(dplyr)
library(usethis)
library(stringr)
library(GenomeInfoDb)

#chr_info<-read.table(textConnection(
#  readLines(
#    gzcon(
#      url('http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/chromInfo.txt.gz')
#    )
#  ))) %>% select(chr=1,length=2) %>% filter(str_length(chr)<6)

#chr_len<-chr_info$length
#names(chr_len)<-chr_info$chr

cytobands_data<-function(genome) {
	chr_info<-Seqinfo(genome = genome)
	chr_info<-chr_info[seqnames(chr_info)[(str_length(seqnames(chr_info)) < 6)]]
	con<-readLines(gzcon( url(str_c('http://hgdownload.cse.ucsc.edu/goldenPath/',genome,'/database/cytoBand.txt.gz'))))
	cytobands<-read.table(textConnection(con),sep='\t',stringsAsFactors = FALSE) %>% select(chr=1,start=2,end=3,name=4,gieStain=5) %>% filter(str_length(chr) < 6) 
	cytobands$start<-cytobands$start+1
	cy<<-cytobands
	ci<<-chr_info
	GRanges(cytobands,seqinfo = chr_info)
}

cytobands_hg18<-cytobands_data('hg18')
usethis::use_data(cytobands_hg18,overwrite=TRUE)
cytobands_hg19<-cytobands_data('hg19')
usethis::use_data(cytobands_hg19,overwrite=TRUE)
cytobands_hg38<-cytobands_data('hg38')
usethis::use_data(cytobands_hg38,overwrite=TRUE)
cytobands_mm10<-cytobands_data('mm10')
usethis::use_data(cytobands_mm10,overwrite=TRUE)
cytobands_mm9<-cytobands_data('mm9')
usethis::use_data(cytobands_mm9,overwrite=TRUE)
cytobands_mm8<-cytobands_data('mm8')
usethis::use_data(cytobands_mm8,overwrite=TRUE)
