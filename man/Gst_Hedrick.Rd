\name{Gst_Hedrick}
\alias{Gst_Hedrick}
\title{Calculate Nei's Gst using estimators for Hs and Ht}
\usage{
  Gst_Hedrick(x)
}
\arguments{
  \item{x}{genind object (from package adegenet)}
}
\description{
  This function calculates Hedrick's G'st from a genind
  object
}
\details{
  Takes a genind object with population information and
  calculates Hedrick's G''st. This Returns a list with
  values for each locus as well as a global estimates

  Because estimators of Hs and Ht are used, it's possible
  to have negative estimates of Gst. You should treat such
  results as zeros (or estimating a value close to zero,
  and getting it a little wrong)
}
\examples{
data(nancycats)
Gst_Hedrick(nancycats)
}
\references{
  Hedrick, PW. (2005), A Standardized Genetic
  Differentiation Measure. Evolution 59: 1633-1638.

  Meirmans PG, Hedrick PW (2011), Assessing population
  structure: FST and related measures. Molecular Ecology
  Resources, 11:5-18
}
\seealso{
  Other diffstat: \code{\link{diff_stats}},
  \code{\link{D_Jost}}, \code{\link{Gst_Nei}}

  Other Hedrick: \code{\link{pairwise_Gst_Hedrick}}
}

