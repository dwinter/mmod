\name{D_Jost}
\alias{D_Jost}
\title{Calculate Jost's D}
\usage{
  D_Jost(x)
}
\arguments{
  \item{x}{genind object (from package adegenet)}
}
\value{
  per.locus values for each D for each locus in the dataset

  global estimtes for D based on overall heterozygosity or
  the harmonic mean of values for each locus
}
\description{
  This function calculates Jost's D from a genind object
}
\details{
  Takes a genind object with population information and
  calculates Jost's D Returns a list with values for each
  locus as well as two global estimates. 'global.het' uses
  the averages of Hs and Ht across all loci while
  'global.harm_mean' takes the harmonic mean of all loci.

  Because estimators of Hs and Ht are used, its possible to
  have negative estimates of D. You should treat these as
  numbers close to zero.
}
\examples{
data(nancycats)
D_Jost(nancycats)
}
\references{
  Jost, L. (2008), GST and its relatives do not measure
  differentiation. Molecular Ecology, 17: 4015-4026.
}
\seealso{
  Other D: \code{\link{pairwise_D}}

  Other diffstat: \code{\link{diff_stats}},
  \code{\link{Gst_Hedrick}}, \code{\link{Gst_Nei}},
  \code{\link{Phi_st_Meirmans}}
}

