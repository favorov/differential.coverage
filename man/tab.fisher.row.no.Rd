% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tabulated.fisher.R
\name{tab.fisher.row.no}
\alias{tab.fisher.row.no}
\title{tab.fisher.row.no}
\usage{
tab.fisher.row.no(Y, N, MY, MN)
}
\arguments{
\item{Y}{sum of column 1 (cases#)}

\item{N}{sum of column 2 (controls#)}

\item{MY}{methylated cases #}

\item{MN}{methylated controls #}
}
\value{
the index of the row
}
\description{
calculates the index of the row (each row is a Fisher table)
that represents the Fisher's table among the matrix that tabulate the Fisher table results
}
