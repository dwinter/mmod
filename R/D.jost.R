D.jost <- function(x){
  n <- length(unique(pop(x)))
  harmN <- harmonic.mean(table(pop(x)))
  D.per.locus <- function(g) {
    #what we need to calculate these stats
    a <- makefreq(genind2genpop(g, quiet=T), quiet=T)[[1]]
    HpS <- sum(1 - apply(a^2, 1, sum, na.rm=TRUE)) / n
    Hs_est <- (2*harmN/(2*harmN-1))*HpS
    HpT <- 1 - sum(apply(a,2,mean, na.rm=TRUE)^2)
    Ht_est <- HpT + Hs_est/(2*harmN*n)
    #The stat itself
    D <- (Ht_est-Hs_est)/(1-Hs_est) * (n/(n-1))
    return(c(Ht_est, Hs_est, D))
  }
 loci <- t(sapply(seploc(x), D.per.locus))
  global_Hs <- mean(loci[,1], na.rm=T)
  global_Ht <- mean(loci[,2], na.rm=T)
  global_D <-  (global_Ht - global_Hs)/(1 - global_Hs ) * (n/(n-1))
  harm_D <- harmonic.mean(loci)
  return(list("per locus"=loci[,3],
              "global (heterozygosities)"=global_D,
              "global (harmonic mean)" = harm_D
        ))

}

