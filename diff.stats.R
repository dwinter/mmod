

diff.stats <- function(x){
  n <- length(unique(pop(x)))
  harmN <- 1/mean(1/table(pop(x)))
  per.locus <- function(g) {
    #what we need to calculate these stats
    a <- makefreq(genind2genpop(g, quiet=T), quiet=T)[[1]]
    HpS <- sum(1 - apply(a^2, 1, sum, na.rm=TRUE)) / n
    Hs_est <- (2*harmN/(2*harmN-1))*HpS
    HpT <- 1 - sum(apply(a,2,mean, na.rm=TRUE)^2)
    Ht_est <- HpT + Hs_est/(2*harmN*n)
    #The stats themselves
    G_est <- (Ht_est-Hs_est)/Ht_est
    D <- (Ht_est-Hs_est)/(1-Hs_est) * (n/(n-1))
    Gprime_st <- G_est * (n-1+Hs_est)/((n-1)*(1-Hs_est))
    #And the results formated as list
    result <- c("Hs" = Hs_est, 
                "Ht" = Ht_est, 
                "Gst_est"=G_est, 
                "G'st" = Gprime_st,
                "D_est" = D)
    return(result)
  }
 loci <- t(sapply(seploc(x), per.locus))
  global_Hs <- mean(loci[,1], na.rm=T)
  global_Ht <- mean(loci[,2], na.rm=T)
  global_G_est <- (global_Ht - global_Hs)/global_Ht
  return(list(loci=loci,
              global=c(Hs =        global_Hs, 
                       Ht =        global_Ht, 
                       Gst_est =   global_G_est, 
                       "G'st_est"= global_G_est*(n-1+global_Hs)/((n-1)*(1-global_Hs)),
                       D_est =     (global_Ht - global_Hs)/(1 - global_Hs ) * (n/(n-1))
                       )
              )
         )
}






#for one locus (allele from makefreq)R
per.locus <- function(g) {
  #what we need to calculate these stats
  n <- length(unique(pop(g)))
  harmN <- 1/mean(1/table(pop(g)))
  a <- makefreq(genind2genpop(g, quiet=T), quiet=T)[[1]]
  HpS <- sum(1 - apply(a^2, 1, sum, na.rm=TRUE)) / n
  Hs_est <- (2*harmN/(2*harmN-1))*HpS
  HpT <- 1 - sum(apply(a,2,mean, na.rm=TRUE)^2)
  Ht_est <- HpT + Hs_est/(2*harmN*n)
  #The stats themselves
  G_est <- (Ht_est-Hs_est)/Ht_est
  D <- (Ht_est-Hs_est)/(1-Hs_est) * (n/(n-1))
  Gprime_st <- G_est * (n-1+Hs_est)/((n-1)*(1-Hs_est))
  #And the results
  result <- c("Hs" = Hs_est, 
              "Ht" = Ht_est, 
              "Gst_est"=G_est, 
              "G'st" = Gprime_st,
              "D_est" = D)
  return(result)
}

genotypic_differntiation <- function(g){
  n <- length(unique(pop(g)))
  harmN <- 1/mean(1/table(pop(g)))
  
  results <- t(sapply(seploc(g), per.locus))
  global_Hs <- mean(results[,1], na.rm=T)
  global_Ht <- mean(results[,2], na.rm=T)
  global_G_est <- (global_Ht - global_Hs)/global_Ht
  global_G_prime_est <- global_G_est*(n-1+global_Hs)/((n-1)*(1-global_Hs))
  global_D <- (global_Ht - global_Hs)/(1 - global_Hs ) * (n/(n-1))
  return(rbind(results,
               c(global_Hs, global_Ht, global_G_est, global_G_prime_est, global_D)  
         ))
}
