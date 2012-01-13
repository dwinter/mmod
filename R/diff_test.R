#' An exact test of population differntiation for Genind objects
#'
#' This function performs and exact test of population differentiation
#' based on allele frequencies in sub-population.  
#'
#' Note, this test returns a p-value _not_ an estimate of effect size.
#' Since most populations have some degree of population differentiation,
#' very large samples are almost guaranteed to return signifcant results.
#' Refer to estimates of D or Gst to ascertain how meaningful such results
#' might be.
#'
#' @param x genind object (from package adegenet)
#' @param sim simulate p-value (required for all but the smallest datasets)
#' @param nreps number of steps used to simulate p-value (default 1000)
#' @export
#' @examples
#' 
#' data(nancycats)
#' diff_test(seploc(nancycats)[[2]], nreps=100)
#' 


diff_test <- function(x, sim=TRUE, nreps=10000){
  # The test to be applied to each locus
  per.locus <- function(locus){
    allele_counts <- locus@tab*2
    if(any(is.na(locus@tab))){
      warning(paste("Ommiting individuals with NAs for locus", locus$loc.names))
      include <- complete.cases(allele_counts)
      allele_counts <- allele_counts[include,]
      pops <- pop(locus)[include, drop=TRUE]
    }
    else {
      pops <- pop(locus)
    }
    # Sum up alleles (rows in matrix) by population using tapply
    pop_alleles <- apply(allele_counts, 2, function(a) tapply(a, pops, sum))
    return(fisher.test(pop_alleles, simulate.p.value=sim, B=nreps)$p.value)
  }

  loci <- (sapply(seploc(x), per.locus))
  return(loci)
}
