#' Calculate differentiation statistics for a genetic dataset
#'
#' This function calculates three different statistics of differentiaion
#' for a genetic dataset. Nei's Gst, Hedrick's Gst and Jost's D
#'
#' @param x genind object (from package adegenet)
#' @export
#' @examples
#' 
#' data(nancycats)
#' diff.stats(nancycats)

diff.stats <- function(x){
  n <- length(unique(pop(x)))
  harmN <- harmonic.mean(table(pop(x)))
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
  return(list("per.locus"=loci,
              global=c(
                Hs = global_Hs, 
                Ht = global_Ht, 
                Gst_est = global_G_est, 
                "G'st_est"= global_G_est*(n-1+global_Hs)/((n-1)*(1-global_Hs)),
                "D_est(het)" = (global_Ht - global_Hs)/(1 - global_Hs ) * (n/(n-1)),
                "D_est(mean)"= harmonic.mean(loci[,5])
              )))
}

