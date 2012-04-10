\name{rgenotypes}
\alias{rgenotypes}
\title{Randomly create genotypes}
\usage{
  rgenotypes(n, ploidy, probs, genind = FALSE,
    pop_name = "A")
}
\arguments{
  \item{n}{integer number of indviduals.}

  \item{ploidy}{integer number of alleles to asign to each
  individual.}

  \item{probs}{vector of probabilies corresponding to
  allele frequences.}

  \item{genind}{boolean if TRUE return a genind object}

  \item{pop_name}{charcter Name for population defined in
  genind object (not required if genind is not TRUE)}
}
\value{
  Either a matrix with individuals in columns, alleles in
  rows or, if genind is TRUE a genind object for one
  population and locus.
}
\description{
  Use the multinomial distribution to randomly create
  genotpes for individuals for given allele frequences. By
  default this function returns a matrix of with alleles in
  rows and individuals in columns. There is an option to
  return a genind object representing the same data (see
  examples).
}
\details{
  Used in \code{\link{chao_bootstrap}}, also exported as it
  may come in handy for other simulations.
}
\examples{
data(nancycats)
obs_allele_freqs <- apply(nancycats$tab[,1:16], 2,mean)
rgenotypes(10, 2, obs_allele_freqs)
obs_allele_freqs_noNA <- apply(nancycats$tab[,1:16], 2,mean, na.rm=TRUE)
rgenotypes(10, 2, obs_allele_freqs_noNA)
}
\seealso{
  \code{\link{rmultinom}} which this function wraps.
}

