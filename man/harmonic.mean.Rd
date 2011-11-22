\name{harmonic.mean}
\alias{harmonic.mean}
\title{Harmonic mean}
\usage{
  harmonic.mean(x)
}
\arguments{
  \item{x}{numeric vector}
}
\description{
  Calculate the harmonic mean of a numeric vector
}
\examples{
data(nancycats)
pop.sizes <- table(pop(nancycats))
harmonic.mean(pop.sizes)
}

