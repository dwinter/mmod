#' Calculate Nei's Gst using estimators for Hs and Ht
#'
#' This function calculates G'st, Hedrick's  correction to Gst 
#' accounting for observed Hs. Nei and Chesser's  estimators of
#' Hs and Ht are used
#'
#' @param x genind object (from package adegenet)
#' @export
#' @examples
#' 
#' data(nancycats) 
#' Gst_Hedrick(nancycats)

Gst_Hedrick <- function(x){
  n <- length(unique(pop(x)))
  harmN <- harmonic_mean(table(pop(x)))
  D.per.locus <- function(g) {
    #what we need to calculate these stats
    a <- makefreq(genind2genpop(g, quiet=T), quiet=T)[[1]]
    HpS <- sum(1 - apply(a^2, 1, sum, na.rm=TRUE)) / n
    Hs_est <- (2*harmN/(2*harmN-1))*HpS
    HpT <- 1 - sum(apply(a,2,mean, na.rm=TRUE)^2)
    Ht_est <- HpT + Hs_est/(2*harmN*n)
    #The stat itself
    G_est <- (Ht_est-Hs_est)/Ht_est
    Gprime_st <- G_est * (n-1+Hs_est)/((n-1)*(1-Hs_est))
    return(c(Hs_est, Ht_est, Gprime_st))
  }
 loci <- t(sapply(seploc(x), D.per.locus))
  global_Hs <- mean(loci[,1], na.rm=T)
  global_Ht <- mean(loci[,2], na.rm=T)
  global_G_est <-  (global_Ht - global_Hs)/global_Ht
  global_Hedrick <-  global_G_est * (n-1+global_Hs)/((n-1)*(1-global_Hs))
  harm_D <- harmonic_mean(loci)
  return(list("per.locus"=loci[,3], "global"=global_Hedrick))

}

