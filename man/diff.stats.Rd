\name{diff.stats}
\alias{diff.stats}
\title{Calculate differentiation statistics for a genetic dataset}
\usage{
  diff.stats(x)
}
\arguments{
  \item{x}{genind object (from package adegenet)}
}
\description{
  This function calculates three different statistics of
  differentiaion for a genetic dataset. Nei's Gst,
  Hedrick's Gst and Jost's D
}
\examples{
data(nancycats)
diff.stats(nancycats)
}

