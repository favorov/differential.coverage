#' differential.coverage: is a collection of routines useful for differential analysis of intervals
#'
#' \code{differential.coverage} is more or less random collection of procedures that are useful for differential analysis of interval sets. 
#' Each set represents a genomwide experiment (probably, MBD, but actually, we just suppose it to be an interval genomewide coverage). 
#' The package provides three classes of functions: prepare probes, test differntial intersection of probes with a collection of intervals, and annotation.
#' 
#' @section Preparing probes:
#'
#' The probes is a set of genomic intervels. It actually can be some pre-given set of genomic intervals (e.g. open chromatin intervals or CpG islands), but if we want to be agnostic, we can prepare a
#' uniform tiling of all the genome. We also can 'inflate' all the probes (intervals).\cr
#' \link{prepare.covering.noodles}\cr
#' \link{inflate.noodles}\cr
#'	
#' @section Test differential intersection of the prbes with a set of intervals:
#' 
#' The probes we have is a common gemtry for all the samples we work with. Also, for each sample, we have set of genomic intervals (lets' call the peaks), e.g. MACS results for a MBDseq,
#' and the samples belong to (two) conditions and for each probe we want to test null hypethesis that the probe interval is equally prone to be intersected with the peaks.\cr
#' \link{count.coverage.of.noodles}\cr
#' \link{indicate.any.coverage.of.noodles}\cr
#' \link{max.peak.score.for.each.noodle}\cr
#' \link{prepare.tabulated.fisher}\cr
#' \link{tab.fisher.row.no}\cr
#' 
#' @section Annotation functions:
#'
#' After we run the tests and filter the probes to have the differentialy covered oned, or for any other interval annotation, we can run the annotation functions to idetify the genes around.
#' The package knows the gene annotatinos of the following genomes: ucsc mm9 and mm10, gencode mm23, mm24 and mm25, ucsc hg 18, hg19, hg38, and gencode hs19, hs26, hs28, hs29, hs32 and hs34.
#' The functions are:\cr
#' \link{closest.gene.by.interval}\cr
#' \link{closest.gene.start.by.interval}\cr
#' \link{genes.intersected}\cr
#' \link{genes.intersected.by.interval}\cr
#' \link{genes.with.TSS.covered}\cr
#' \link{genes.with.TSS.covered.by.interval}\cr
#' \link{get.Known.Gene.List}\cr
#' \link{preceded.gene.by.interval}\cr


#' @docType package
#' @name differential.coverage 
NULL
