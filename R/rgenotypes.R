#' Randomly create genotypes
#'
#' Use the multinomial distribution to randomly create genotpes for individuals
#' for given allele frequences. 
#' 
#' Used in \code{\link{chao_bootstrap}}, also exported as it may come in handy
#' for other simulations
#' 
#' @param n integer number of indviduals 
#' @param ploidy integer number of alleles to asign to each indivudal
#' @param probs vector of probabilies corresponding to allele frequences
#'
#' @return A matrix with individuals in columns, alleles in rows
#' @seealso \code{\link{rmultinom}} which this function wraps
#' @export
#' @examples
#' 
#' data(nancycats)
#' obs_allele_freqs <- apply(nancycats$tab[,1:16], 2,mean)
#' rgenotypes(10, 2, obs_allele_freqs)
#' obs_allele_freqs_noNA <- apply(nancycats$tab[,1:16], 2,mean, na.rm=TRUE)
#' rgenotypes(10, 2, obs_allele_freqs_noNA)

rgenotypes <- function(n, ploidy, probs){
 if(all(is.na(probs))){ 
  return(matrix(rep(probs, n), ncol=n)) 
  }
 else
 return( rmultinom(n, ploidy, probs) ) 
}

