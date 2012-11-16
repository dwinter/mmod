\name{as.genind.DNAbin}
\alias{as.genind.DNAbin}
\title{as.genind.DNAbin}
\usage{
  as.genind.DNAbin(x, pops)
}
\arguments{
  \item{x}{object of class DNAbin}

  \item{pops}{vector of population assignemnts for each
  sequence}
}
\value{
  genind
}
\description{
  Convert a DNAbin object into a genind object
}
\examples{
library(pegas)
data(woodmouse)
wm <- as.genind.DNAbin(woodmouse, rep(c("A", "B", "C"), each=5))
diff_stats(wm)
}

