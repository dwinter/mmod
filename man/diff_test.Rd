\name{diff_test}
\alias{diff_test}
\title{An exact test of population differntiation for Genind objects}
\usage{
  diff_test(x, sim = TRUE, nreps = 10000)
}
\arguments{
  \item{x}{genind object (from package adegenet)}
}
\description{
  This function performs and exact test of population
  differentiation based on allele frequencies in
  sub-population.
}
\details{
  Note, this test returns a p-value _not_ an estimate of
  effect size. Since most populations have some degree of
  population differentiation, very large samples are almost
  guaranteed to return signifcant results. Refer to
  estimates of D or Gst to ascertain how meaningful such
  results might be.
}
\examples{
data(nancycats)
diff_test(nancycats, nreps=100)
}

