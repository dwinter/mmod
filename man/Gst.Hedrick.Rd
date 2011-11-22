\name{Gst.Hedrick}
\alias{Gst.Hedrick}
\title{Calculate Nei's Gst using estimators for Hs and Ht}
\usage{
  Gst.Hedrick(x)
}
\arguments{
  \item{x}{genind object (from package adegenet)}
}
\description{
  This function calculates G'st, Hedrick's correction to
  Gst accounting for observed Hs. Nei and Chesser's
  estimators of Hs and Ht are used
}
\examples{
data(nancycats)
Gst.Hedrick(nancycats)
}

