\name{syntable_db}
\alias{syntable_db}
\title{Synoptic tables for long-format vegetation databases}
\usage{
syntable_db(db, cluster, abund = "percentage", type = "percfreq", digits = 0)
}
\arguments{
  \item{db}{Data frame in long format containing at least the columns \code{RELEVE_NR}, \code{TaxonName} and \code{Cover_Perc}.}
  \item{cluster}{Integer or character vector/factor with classification cluster identity for each \code{RELEVE_NR}.}
  \item{abund}{Type of abundances. Define whether \code{Cover_Perc} contains percentage cover (\code{abund = "percentage"}, default) or presence/absence data (\code{abund = "pa"}, with values 0/1).}
  \item{type}{Type of synoptic table output \code{type = c("percfreq", "totalfreq", "mean", "median", "phi")}.}
  \item{digits}{Integer indicating the number of decimal places to be displayed in result tables (default 0).}
}
\value{
An (invisible) list with components:
\item{\$syntable}{unordered synoptic table for given species and clusters}
\item{\$samplesize}{total number of samples per cluster}
}
\seealso{\code{\link{syntable}}, \code{\link{synsort}}}
\description{
\code{syntable_db} calculates synoptic tables for long-format vegetation data where species occur in rows and samples are given by \code{RELEVE_NR}.
}
