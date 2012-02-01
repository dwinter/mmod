\name{chao_bootstrap}
\alias{chao_bootstrap}
\title{Produce bootstrap samples from each subpopulation of a genind object}
\usage{
  chao_bootstrap(x, nreps = 1000)
}
\arguments{
  \item{x}{genind object (from package adegenet)}

  \item{nreps}{numeric number of bootstrap replicates to
  perform (default 1000)}
}
\description{
  This function produces bootstrap samples from a genind
  object, with each subpopulation resampled according to
  its size. Because there are many statistics that you may
  wish to calculte from these samples, this function
  returns a list of genind objects representing bootsrap
  samples that can then be futher processed (see examples).
}
\examples{
\dontrun{
data(nancycats)
bs <- chao_bootstrap(nancycats)
bs_D <- summarise_bootstrap(bs, D_Jost)
}
}
\seealso{
  Other resample: \code{\link{jacknife_populations}},
  \code{\link{summarise_bootsrap}}
}

