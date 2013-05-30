permute.HsHt <- function(x){
  temp <- x
  pop(temp) <- sample(pop(x))
  n <- length(unique(pop(x)))
  harmN <- harmonic.mean(table(pop(x)))
  H.per.locus <- function(g) {
    #what we need to calculate these stats
    a <- makefreq(genind2genpop(g, quiet=T), quiet=T)[[1]]
    HpS <- sum(1 - apply(a^2, 1, sum, na.rm=TRUE)) / n
    Hs_est <- (2*harmN/(2*harmN-1))*HpS
    HpT <- 1 - sum(apply(a,2,mean, na.rm=TRUE)^2)
    Ht_est <- HpT + Hs_est/(2*harmN*n)
    #The stat itself
    return(c(Ht = Ht_est, Hs = Hs_est))
  }
 loci <- t(sapply(seploc(temp), H.per.locus))
  return(list("per locus"=loci,
              "global" =c(Ht = mean(loci[,1], na.rm=T),
                          Hs = mean(loci[,2], na.rm=T)
             
         )))

}

