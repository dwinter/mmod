\name{diff_stats}
\alias{diff_stats}
\title{Calculate differentiation statistics for a genind objects}
\usage{
  diff_stats(x)
}
\arguments{
  \item{x}{genind object (from package adegenet)}
}
\description{
  This function calculates three different statistics of
  differentiation for a genetic dataset. Nei's Gst,
  Hedrick's G'st and Jost's D
}
\details{
  See individual functions D_Jost(), Gst_Hedrick() and
  Gst_Nei() for more details
}
\examples{
data(nancycats)
diff_stats(nancycats)
}
\references{
  Hedrick, PW. (2005), A Standardized Genetic
  Differentiation Measure. Evolution 59: 1633-1638.

  Jost, L. (2008), GST and its relatives do not measure
  differentiation. Molecular Ecology, 17: 4015-4026.

  Meirmans PG, Hedrick PW (2011), Assessing population
  structure: FST and related measures. Molecular Ecology
  Resources, 11:5-18

  Nei M. (1973) Analysis of gene diversity in subdivided
  populations. PNAS: 3321-3323.

  Nei M, Chesser RK. (1983). Estimation of fixation indices
  and gene diversities. Annals of Human Genetics. 47:
  253-259.
}
\seealso{
  Other diffstat: \code{\link{D_Jost}},
  \code{\link{Gst_Hedrick}}, \code{\link{Gst_Nei}}
}

