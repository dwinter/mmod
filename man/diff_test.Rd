\name{diff_test}
\alias{diff_test}
\title{An exact test of population differntiation for Genind objects}
\usage{
  diff_test(x, sim = TRUE, nreps = 10000)
}
\arguments{
  \item{x}{genind object (from package adegenet)}

  \item{sim}{simulate p-value (required for all but the
  smallest datasets)}

  \item{nreps}{number of steps used to simulate p-value
  (default 1000)}
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
diff_test(seploc(nancycats)[[2]], nreps=100)
}

