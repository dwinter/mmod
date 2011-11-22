\name{jacknife.populations}
\alias{jacknife.populations}
\title{Calculate differentiation stats for a jacknife sample of a Genind opject}
\usage{
  jacknife.populations(x, sample_frac = 0.5, nreps = 1000)
}
\arguments{
  \item{x}{genind object (from package adegenet)}

  \item{sample_frac}{fraction of pops to sample in each
  replication (default 0.5)}

  \item{nreps}{number of jacknife replicates to run
  (default 1000)}
}
\description{
  Calcutes
}
\examples{
\dontrun{
data(nancycats)
jacknife.populations(nancycats)
}
}

