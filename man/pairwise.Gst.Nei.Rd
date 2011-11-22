\name{pairwise.Gst.Nei}
\alias{pairwise.Gst.Nei}
\title{Calculates pairwise values of Nei's Gst}
\usage{
  pairwise.Gst.Nei(x)
}
\arguments{
  \item{x}{genind object (from package adegenet)}
}
\description{
  This function calculates Nei's Gst, a measure of genetic
  differentiation, between all combinations of populaitons
  in a genind object.
}
\examples{
data(nancycats)
pairwise.Gst.Nei(nancycats[1:26,])
}

