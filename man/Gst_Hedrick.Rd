\name{Gst_Hedrick}
\alias{Gst_Hedrick}
\title{Calculate Nei's Gst using estimators for Hs and Ht}
\usage{
  Gst_Hedrick(x)
}
\arguments{
  \item{x}{genind object (from package adegenet)}
}
\description{
  This function calculates Hedrick's G'st from a genind
  object
}
\details{
  Takes a genind object with population information and
  calculates Hedrick's G'st. This Returns a list with
  values for each locus as well as a global estimates

  Because estimators of Hs and Ht are used, its possible to
  have negative estimates of Gst. You should treat such
  results as zeros (or estimating a value close to zero,
  and getting it a little wrong)
}
\examples{
data(nancycats)
Gst_Hedrick(nancycats)
}

