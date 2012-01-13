\name{pairwise_D}
\alias{pairwise_D}
\title{Calculates pairwise values of Jost's D}
\usage{
  pairwise_D(x)
}
\arguments{
  \item{x}{genind object (from package adegenet)}
}
\description{
  This function calculates Jost's D, a measure of genetic
  differentiation, between all combinations of populaitons
  in a genind object.
}
\examples{
data(nancycats)
pairwise_D(nancycats[1:26,])
}

