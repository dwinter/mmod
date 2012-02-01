\name{rgenotypes}
\alias{rgenotypes}
\title{Randomly create genotypes}
\usage{
  rgenotypes(n, ploidy, probs)
}
\arguments{
  \item{n}{integer number of indviduals}

  \item{ploidy}{integer number of alleles to asign to each
  indivudal}

  \item{probs}{vector of probabilies corresponding to
  allele frequences}
}
\value{
  A matrix with individuals in columns, alleles in rows
}
\description{
  Use the multinomial distribution to randomly create
  genotpes for individuals for given allele frequences.
}
\details{
  Used in \code{\link{chao_bootstrap}}, also exported as it
  may come in handy for other simulations
}
\examples{
data(nancycats)
obs_allele_freqs <- apply(nancycats$tab[,1:16], 2,mean)
rgenotypes(10, 2, obs_allele_freqs)
obs_allele_freqs_noNA <- apply(nancycats$tab[,1:16], 2,mean, na.rm=TRUE)
rgenotypes(10, 2, obs_allele_freqs_noNA)
}
\seealso{
  \code{\link{rmultinom}} which this function wraps
}

