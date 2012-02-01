#' Calculate Jost's D 
#'
#' This function calculates Jost's D from a genind object
#'
#' Takes a genind object with population information and calculates Jost's D
#' Returns a list with values for each locus as well as two global estimates.
#' 'global.het' uses the averages of Hs and Ht across all loci while 
#' 'global.harm_mean' takes the harmonic mean of all loci.
#' 
#' Because estimators of Hs and Ht are used, its possible to have negative
#' estimates of D. You should treat these as numbers close to zero.
#'
#' @param x genind object (from package adegenet)
#' @export
#' @examples
#' 
#' data(nancycats)
#' D_Jost(nancycats)
#' @references
#'  Jost, L. (2008), GST and its relatives do not measure differentiation. Molecular Ecology, 17: 4015-4026. 
#' @family diffstat
#' @family D

D_Jost <- function(x){
  n <- length(unique(pop(x)))
  harmN <- harmonic_mean(table(pop(x)))
  pops <- pop(x)
  D.per.locus <- function(g) {
    a <- apply(g@tab,2,function(row) tapply(row, pops, mean, na.rm=TRUE))
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
  harm_D <- harmonic_mean(loci[,3])
  return(list("per.locus"=loci[,3],
              "global.het"=global_D,
              "global.harm_mean" = harm_D
        ))

}

