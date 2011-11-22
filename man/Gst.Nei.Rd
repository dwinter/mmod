\name{Gst.Nei}
\alias{Gst.Nei}
\title{Calculate Nei's Gst using estimators for Hs and Ht}
\usage{
  Gst.Nei(x)
}
\arguments{
  \item{x}{genind object (from package adegenet)}
}
\description{
  This function calculates Gst following Nei's method and
  using Nei and Chesser's estimators for Hs and Ht
}
\examples{
data(nancycats)
Gst.Nei(nancycats)
}

