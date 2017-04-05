#' Apply a differentiation statistic to a bootstrap sample
#'
#' This function applies a differentiation statistic (eg, D_Jost, Gst_Hedrick or 
#' Gst_Nei) to a list of genind objects, possibly produced with
#' chao_bootsrap or jacknife_populations. 
#' 
#' Two different approaches are used for calculating confidence intervals in the
#' results. The estimates given by \code{lower.percentile} and \code{upper.percentile}
#' are simply the \code{2.5}th and \code{97.5}th precentile of the statistic
#' across bootstrap samples. Note, the presence or rare alleles in some
#' populations can bias bootstrapping procedures such that these intervals
#' are not centered on the observed value. The mean of statistic across
#' samples is returned as \code{mean.bs} and can be used to correct biased
#' bootsrap samples. Alternatively, \code{lower.normal} and \code{upper.normal}
#' form a confidence interval centered on the observed value of the statistic
#' and using the standard deviation of the statistic across replicates to
#' generate limits (sometimes called the normal-method of obtaining a confidence
#' interval). The print function for objects returned by this function displays 
#' the normal-method confidence intervals.
#' 
#'
#' @param bs list of genind objects
#' @param statistic differentiation statistic to apply (the function itself, 
#' as with apply family functions)
#' @family resample
#' @return per.locus:  \code{matirx} of statistics calculated for each locus (column) and each 
#' bootstrap replicate (row).
#' @return global.het: \code{vector} of global estimates calculated from overall 
#' heterozygosity 
#' @return global.het: \code{vector} of global estimates calculated from harmonic
#' mean of statistic (only applied to D_Jost)
#' @return summary.loci: \code{data.frame} summarising the distribution of the
#' chosen statistic across replicates. Details of the different confidence
#' intervals are given in details
#' @return summary.global_het: A vector containing the same measures as
#' \code{summary.loci} but for a global value of the statistic calculated from
#' all loci
#' @return summary.global_harm: As with \code{summary.global_het} but calculated
#' from the harmonic mean of the statistic across loci (only applies to D_Jost)
#' @importFrom stats quantile
#' @importFrom stats sd
#' @export
#' @examples
#'\dontrun{  
#' data(nancycats)
#' bs <- chao_bootstrap(nancycats)
#' summarise_bootstrap(bs, D_Jost)
#'}



summarise_bootstrap <- function(bs, statistic){
    obs <- statistic(bs$obs)
    nreps <- length(bs)
    stats <- sapply(bs$BS, statistic)
    loc_stats <- do.call(rbind, stats["per.locus",])

    res <-list("per.locus"= loc_stats,
               "global.het"=unlist(stats[2,])
    )
    nloc <- length(obs$per.locus)
    res$summary.loci <- data.frame(
        "locus"    = names(obs$per.locus),                           
        "observed" = obs$per.locus,
        "lower.normal" = numeric(nloc), 
        "upper.normal" = numeric(nloc), 
        "std.dev" = numeric(nloc), 
        "mean.bs" = numeric(nloc), 
        "lower.percentile" = numeric(nloc), 
        "upper.percentile" = numeric(nloc),
        stringsAsFactors=FALSE #hell no
    )
    for(i in 1:nloc){
        res$summary.loci[i,2:8] <- combined_summary_stats(obs$per.locus[i], bs_summary_stats(loc_stats[,i]))
    }
                                
    if(identical(statistic, D_Jost)){
        res$global.harm <- unlist(stats[3,])
        if(any(is.na(res$global.harm))){
            res$global.harm[is.na(res$global.harm)] <- 0
            warning("Bootstrap distribution of D_Jost includes negative values, harmonic mean is undefined")
        }
        res$summary.global.harm <- combined_summary_stats(obs$global.harm_mean, bs_summary_stats(res$global.harm))
        res$summary.global.het<- combined_summary_stats(obs$global.het, bs_summary_stats(res$global.het))
    } else {
        res$summary.global.het <- combined_summary_stats(obs$global, bs_summary_stats(res$global.het))
    }
    class(res) <- "summarised_bs"
    res
}

#' @export 

print.summarised_bs <- function(x, ...){
  
  print_line <- function(x) sprintf("%.4f\t(%.3f-%.3f)\n", x["observed"], x["lower.normal"], x["upper.normal"]) 

  cat("\nEstimates for each locus\n")
  cat("Locus\tMean\t 95% CI\n")
  for(i in 1:dim(x$summary.loci)[1]){
    cat( x$summary.loci$locus[i], print_line(x$summary.loci[i,]), sep="\t")    
  }
  cat("\nGlobal Estimate based on average heterozygosity\n")
  cat(print_line(x$summary.global.het))
  if(!is.null(x$summary.global.harm)){
    cat("\nGlobal Estimate based on harmonic mean of statistic\n")
    cat(print_line(x$summary.global.harm))
  }
}

bs_summary_stats <- function(B){
    
    structure(c(sd(B), mean(B), quantile(B, c(0.025, 0.975))), 
              .Names=c("std.dev", "mean", "lower.percentile", "upper.precentile"))
}

combined_summary_stats <- function(obs, summ){
    c(observed=obs, lower.normal=obs - summ[["std.dev"]]*1.96, upper.normal=obs + summ[["std.dev"]]*1.96, summ)
}
