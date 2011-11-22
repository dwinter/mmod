D.Jost <- function(x){
  n <- length(unique(pop(x)))
  harmN <- harmonic.mean(table(pop(x)))
  D.per.locus <- function(g) {
    
    a <- makefreq(genind2genpop(g, quiet=T), quiet=T)[[1]]
    HpS <- sum(1 - apply(a^2, 1, sum, na.rm=TRUE)) / n
    Hs_est <- (2*harmN/(2*harmN-1))*HpS
    HpT <- 1 - sum(apply(a,2,mean, na.rm=TRUE)^2)
    Ht_est <- HpT + Hs_est/(2*harmN*n)
    
    D <- (Ht_est-Hs_est)/(1-Hs_est) * (n/(n-1))
    return(c(Hs_est, Ht_est, D))
  }
 loci <- t(sapply(seploc(x), D.per.locus))
  global_Hs <- mean(loci[,1], na.rm=T)
  global_Ht <- mean(loci[,2], na.rm=T)
  global_D <-  (global_Ht - global_Hs)/(1 - global_Hs ) * (n/(n-1))
  harm_D <- harmonic.mean(loci)
  return(list("per.locus"=loci[,3],
              "global.het"=global_D,
              "global.harm_mean" = harm_D
        ))

}

