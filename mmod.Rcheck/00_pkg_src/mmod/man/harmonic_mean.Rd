\name{harmonic_mean}
\alias{harmonic_mean}
\title{Harmonic mean}
\usage{
  harmonic_mean(x)
}
\arguments{
  \item{x}{numeric vector}
}
\description{
  Calculate the harmonic mean of a numeric vector (will
  return NA if there are any negative numbers in the
  vector)
}
\examples{
data(nancycats)
pop.sizes <- table(pop(nancycats))
harmonic_mean(pop.sizes)
}

