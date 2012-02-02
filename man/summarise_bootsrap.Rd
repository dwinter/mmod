\name{summarise_bootsrap}
\alias{summarise_bootsrap}
\title{Apply a differentiation statistic to a bootstrap sample}
\usage{
  summarise_bootsrap(bs, statistic)
}
\arguments{
  \item{bs}{list of genind objects}

  \item{statistic}{differentiation statistic to apply.}
}
\description{
  This function applies a differentiation statistic (eg,
  D_Jost, Gst_Hedrick or Gst_Nei) to a list of genind
  objects, possibly produced with chao_bootsrap or
  jacknife_populations. The resulting list contains a
  matrix of values with the statistic for each locus as
  well as a global estimate for every object in the sample.
  Additionally, mean and 95% confidence intervals are
  calculated for each set of statisics A custom print
  method that displays these summaries is provided.
}
\examples{
\dontrun{
data(nancycats)
bs <- chao_bootstrap(nancycats)
summarise_bootstrap(bs, D_Jost)
}
}
\seealso{
  Other resample: \code{\link{chao_bootstrap}},
  \code{\link{jacknife_populations}}
}

