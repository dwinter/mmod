#' Calculate differentiation statistics for a genind objects
#'
#' This function calculates three different statistics of differentiation
#' for a genetic dataset. Nei's Gst, Hedrick's G'st and Jost's D
#' 
#' See individual functions D_Jost(), Gst_Hedrick() and Gst_Nei() for more
#' details
#'
#' @param x genind object (from package adegenet)
#' @export
#' @examples
#' 
#' data(nancycats)
#' diff_stats(nancycats)
#' @references
#'  Hedrick, PW. (2005), A Standardized Genetic Differentiation Measure. Evolution 59: 1633-1638. 
#' @references
#'  Jost, L. (2008), GST and its relatives do not measure differentiation. Molecular Ecology, 17: 4015-4026.
#' @references
#'  Nei M. (1973) Analysis of gene diversity in subdivided populations. PNAS: 3321-3323. 
#' @references
#'  Nei M, Chesser RK. (1983). Estimation of fixation indices and gene diversities. Annals of Human Genetics. 47: 253-259.
#' @family diffstat


diff_stats <- function(x){
  n <- length(unique(pop(x)))
  harmN <- harmonic_mean(table(pop(x)))
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
                "Gst"=G_est, 
                "Gprime_st" = Gprime_st,
                "D" = D)
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
                "Gprime_st"= global_G_est*(n-1+global_Hs)/((n-1)*(1-global_Hs)),
                "D_het" = (global_Ht - global_Hs)/(1 - global_Hs ) * (n/(n-1)),
                "D_mean"= harmonic_mean(loci[,5])
              )))
}

