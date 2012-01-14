\name{Gst_Nei}
\alias{Gst_Nei}
\title{Calculate Nei's Gst using estimators for Hs and Ht}
\usage{
  Gst_Nei(x)
}
\arguments{
  \item{x}{genind object (from package adegenet)}
}
\description{
  This function calculates Gst following Nei's method and
  using Nei and Chesser's estimators for Hs and Ht
}
\examples{
data(nancycats)
Gst_Nei(nancycats)
}
\references{
  Nei M. (1973) Analysis of gene diversity in subdivided
  populations. PNAS: 3321-3323.

  Nei M, Chesser RK. (1983). Estimation of fixation indices
  and gene diversities. Annals of Human Genetics. 47:
  253-259.
}
\seealso{
  Other diffstat: \code{\link{diff_stats}},
  \code{\link{D_Jost}}, \code{\link{Gst_Hedrick}}

  Other Nei: \code{\link{pairwise_Gst_Nei}}
}

