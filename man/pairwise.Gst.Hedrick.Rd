\name{pairwise.Gst.Hedrick}
\alias{pairwise.Gst.Hedrick}
\title{Calculates pairwise values of Hedrick's G'st}
\usage{
  pairwise.Gst.Hedrick(x)
}
\arguments{
  \item{x}{genind object (from package adegenet)}
}
\description{
  This function calculates Hedrick's G'st, a measure of
  genetic differentiation, between all combinations of
  populaitons in a genind object.
}
\examples{
data(nancycats)
pairwise.Gst.Hedrick(nancycats[1:26,])
}

