#' An exact test of population differntiation for Genind objects
#'
#' This function performs and exact test of population differentiation
#' based on allele frequencies between population. 
#'
#' See individual functions D_Jost(), Gst_Hedrick() and Gst_Nei() for more
#' details
#'
#' @param x genind object (from package adegenet)
#' @export
#' @examples
#' 
#' data(nancycats)
#' diff_stats(nancycats)
#' @references
#'  Jost, L. (2008), GST and its relatives do not measure differentiation. Molecular Ecology, 17: 4015-4026. 




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

  loci <- t(sapply(seploc(x), per.locus))
  return(loci)
}
