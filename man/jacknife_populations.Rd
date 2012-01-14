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
  from a Genind object and calculates differentiation stats
  for each sample.
}
\examples{
\dontrun{
data(nancycats)
obs <- diff_stats(nancycats)
jn <- jacknife_populations(nancycats)
D_sampled <-jn[5,]
hist(D_sampled)
abline(h=obs$global)
}
}

