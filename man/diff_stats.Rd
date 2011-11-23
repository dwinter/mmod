\name{diff_stats}
\alias{diff_stats}
\title{Calculate differentiation statistics for a genind objects}
\usage{
  diff_stats(x)
}
\arguments{
  \item{x}{genind object (from package adegenet)}
}
\description{
  This function calculates three different statistics of
  differentiaion for a genetic dataset. Nei's Gst,
  Hedrick's G'st and Jost's D
}
\examples{
data(nancycats)
diff_stats(nancycats)
}
\references{
  Jost, L. (2008), GST and its relatives do not measure
  differentiation. Molecular Ecology, 17: 4015–4026
}

