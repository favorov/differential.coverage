% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/noodles.coverage.R
\name{max.peak.score.for.each.noodle}
\alias{max.peak.score.for.each.noodle}
\title{max.peak.score.for.each.noodle}
\usage{
\method{max}{peak.score.for.each.noodle}(noodles, bedfilnames, bed.ids = bedfilnames)
}
\arguments{
\item{noodles}{\code{GRanges} with the intervals}

\item{bedfilnames}{list of names of bedfiles, one per sample, with some (e.g. methylation) coverage information}

\item{bed.ids}{optional list of names for the samples, they will be used as column names in the result. Default = \code{bedfilnames}}
}
\value{
\code{Matrix}, each row correspond to a noodle; columns are samples, sparse=TRUE, value=max(score) or 0
}
\description{
Could be an alterantive central function of all the differential.coverage package. It gets a set of intervals (noodles) as a \code{GRanges} and a list of names of bed files.
Each bedfile represents a sample. For each noodle and each sample, the maximal score of a peak that intersect this noodle in the sample, if any, or 0 othewise, is returned. 
Supposed to called for annotation of a few selected noodles.
}
\seealso{
\code{differential.coverage}
}
