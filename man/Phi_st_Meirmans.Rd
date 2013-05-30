\name{Phi_st_Meirmans}
\alias{Phi_st_Meirmans}
\title{Calculate Phi_st from a genind object}
\usage{
  Phi_st_Meirmans(x)
}
\arguments{
  \item{x}{genind object (from package adegenet)}
}
\value{
  per.locus Phi_st estimate for each locus

  global Phi_st estimate across all loci
}
\description{
  This function calculates Meirmans' corrected version of
  Phi_st, an Fst analog produced using the AMOVA framework.
  Note, the global estimate produced by this function is
  calculated as the mean distance between individuals
  across all loci, and this exlcuded individuals with one
  or more missing value.
}
\examples{
data(nancycats)
Phi_st_Meirmans(nancycats[1:26,])
}
\references{
  Meirmans, PW. (2005), Using the AMOVA framework to
  estimate a standardized genetic differentiation measure.
  Evolution 60: 2399-402.

  Excoffier, L., Smouse, P., Quattro, J. (1992), Analysis
  of molecular variance inferred from metric distances
  among DNA haplotypes: application to human mitochondrial
  DNA restriction data. Genetics 131: 479-91
}
\seealso{
  Other diffstat: \code{\link{D_Jost}},
  \code{\link{Gst_Hedrick}}, \code{\link{Gst_Nei}},
  \code{\link{diff_stats}}
}

