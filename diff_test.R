diff_test <- function(x, test="exact", sim=TRUE, nreps=10000){

  per.locus <- function(locus){
    allele_counts <- locus@tab*2
    if(any(is.na(locus@tab))){
      warning("Ommiting individuals with NAs")
      include <- complete.cases(allele_counts)
      allele_counts <- allele_counts[include,]
      pops <- pop(locus)[include]
    }
    else {
      pops <- pop(locus)
    }
    
    pop_alleles <- apply(allele_counts, 2, function(a) tapply(a, pops, sum))
    return(fisher.test(pop_alleles, simulate.p.value=sim, B=nreps)$p.value)
  }

  loci <- t(sapply(seploc(x), per.locus))
  return(loci)
}
