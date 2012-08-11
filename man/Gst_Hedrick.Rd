\name{Gst_Hedrick}
\alias{Gst_Hedrick}
\title{Calculate Nei's Gst using estimators for Hs and Ht}
\usage{
  Gst_Hedrick(x)
}
\arguments{
  \item{x}{genind object (from package adegenet)}
}
\value{
  per.locus values for each G''st for each locus in the
  dataset

  global estimtes for G''st based on overall heterozygosity
}
\description{
  This function calculates Hedrick's G'st from a genind
  object
}
\details{
  Takes a genind object with population information and
  calculates Hedrick's G''st.

  Because estimators of Hs and Ht are used, it's possible
  to have negative estimates of G''st. You should treat
  such results as zeros (or an attempt to estimate a very
  low number with some error which might push it below
  zero)
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
  Other Hedrick: \code{\link{pairwise_Gst_Hedrick}}

  Other diffstat: \code{\link{D_Jost}},
  \code{\link{Gst_Nei}}, \code{\link{Phi_st_Meirmans}},
  \code{\link{diff_stats}}
}

