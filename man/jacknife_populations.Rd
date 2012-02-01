\name{jacknife_populations}
\alias{jacknife_populations}
\title{Calculate differentiation stats for a jacknife sample of a Genind opject}
\usage{
  jacknife_populations(x, sample_frac = 0.5, nreps = 1000)
}
\arguments{
  \item{x}{genind object (from package adegenet)}

  \item{sample_frac}{fraction of pops to sample in each
  replication (default 0.5)}

  \item{nreps}{number of jacknife replicates to run
  (default 1000)}
}
\description{
  Makes a series of jacknife samples across populations
  from a Genind object. This functino resurns a list of
  genind objects that can then be further processed (see
  examples).
}
\examples{
\dontrun{
data(nancycats)
obs <- diff_stats(nancycats)
jn <- jacknife_populations(nancycats)
jn.D <- summarise_bootstrap(jn, D_Jost)
}
}
\seealso{
  Other resample: \code{\link{chao_bootstrap}},
  \code{\link{summarise_bootsrap}}
}

