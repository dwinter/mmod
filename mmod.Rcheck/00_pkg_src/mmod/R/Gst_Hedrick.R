#' Calculate Nei's Gst using estimators for Hs and Ht
#'
#' This function calculates Hedrick's G'st from a genind object
#'
#' Takes a genind object with population information and calculates Hedrick's 
#' G'st. This Returns a list with values for each locus as well as a global estimates
#' 
#' Because estimators of Hs and Ht are used, it's possible to have negative
#' estimates of Gst. You should treat such results as zeros (or estimating a
#' value close to zero, and getting it a little wrong)
#'
#' @param x genind object (from package adegenet)
#' @export
#' @examples
#' 
#' data(nancycats) 
#' Gst_Hedrick(nancycats)
#' @references
#'  Hedrick, PW. (2005), A Standardized Genetic Differentiation Measure. Evolution 59: 1633-1638. 
#' @family diffstat
#' @family Hedrick

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

