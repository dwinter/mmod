#' Calculate Nei's Gst using estimators for Hs and Ht
#'
#' This function calculates Gst following Nei's method and using
#' Nei and Chesser's estimators for Hs and Ht
#'
#' @param x genind object (from package adegenet)
#' @export
#' @examples
#' 
#' data(nancycats)
#' Gst_Nei(nancycats)
#' @references
#'  Nei M. (1973) Analysis of gene diversity in subdivided populations. PNAS: 3321-3323. 
#' @references
#'  Nei M, Chesser RK. (1983). Estimation of fixation indices and gene diversities. Annals of Human Genetics. 47: 253-259.
#' @family diffstat
#' @family Nei

Gst_Nei <- function(x){
  n <- length(unique(pop(x)))
  harmN <- harmonic_mean(table(pop(x)))
  Gst.per.locus <- function(g) {
    #what we need to calculate these stats
    a <- makefreq(genind2genpop(g, quiet=T), quiet=T)[[1]]
    HpS <- sum(1 - apply(a^2, 1, sum, na.rm=TRUE)) / n
    Hs_est <- (2*harmN/(2*harmN-1))*HpS
    HpT <- 1 - sum(apply(a,2,mean, na.rm=TRUE)^2)
    Ht_est <- HpT + Hs_est/(2*harmN*n)
    #The stat itself
    G_est <- (Ht_est-Hs_est)/Ht_est
    return(c(Hs_est, Ht_est, G_est))
  }
  
 loci <- t(sapply(seploc(x), Gst.per.locus))
  global_Hs <- mean(loci[,1], na.rm=T)
  global_Ht <- mean(loci[,2], na.rm=T)
  global_G_est <-  (global_Ht - global_Hs)/global_Ht
  
  return(list("per.locus"=loci[,3], "global"=global_G_est))

}

