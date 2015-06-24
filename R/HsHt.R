# Calculate estimators for Hs and Ht
#
# This function calculates Nei and Chesser's estimators for Hs and Ht
# 
# Note, the individuals functions in the diff_stat family (listed below)
# return these estimates (becuase they use this function internally). This
# function is not exported (if someone has a use for it I can change this).

HsHt <- function(x){
    #start by dropping missing for this locus
    x <- x[apply(x@tab, 1, function(i) any(i > 0)),]
    n_by_pop <-  table(pop(x))
    pops <- pop(x)
    n <- nPop(x)
    harmN <- harmonic_mean(n_by_pop[n_by_pop >0])
    a <- apply(apply(x@tab, 2, "/", ploidy(x)), 2, function(L) 
                   tapply(L, pops, mean, na.rm=TRUE))
    HpS <- sum(1 - apply(a^2, 1, sum, na.rm=TRUE)) / n
    Hs_est <- (2*harmN/(2*harmN-1))*HpS
    HpT <- 1 - sum(apply(a,2,mean, na.rm=TRUE)^2)
    Ht_est <- HpT + Hs_est/(2*harmN*n)
    return(c(Ht_est, Hs_est, n))
}
