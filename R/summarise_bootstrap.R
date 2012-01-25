summarise_bootsrap <- function(bs, statistic){
  nreps <- length(bs)
  stats <- sapply(bs, statistic)
  loc_stats <- matrix(unlist(stats[1,]), nrow=nreps, 
               dimnames=list(paste("rep", 1:nreps, sep=""), bs[[1]]@loc.names)
               )
  res <-list("per.locus"= loc_stats,
             "global.het"=unlist(stats[2,])
            )
  #Only D_Jost has another estimate of global value
  if(identical(statistic, D_Jost)){
    res$global.harm <- unlist(stats[3,])
    }
  class(res) <- "summarised_bs"
  return(res)
}

print.summarised_bs <- function(x, digits=4){
  
  summarise <- function(x){
  res <- round(c(mean=mean(x), quantile(x, c(0.025, 0.975))), digits)
  return(paste("\t", res[1], "\t(", res[2],"-", res[3], ")\n", sep="")) 
  }
  
  loc.names <- colnames(x$per.locus)
  loc.results <- apply(x$per.locus, 1, summarise)
  cat("\nEstimates for each locus\n")
  cat("Locus\tMean\t 95% CI\n")
  for(i in 1:length(loc.names)){
    cat(paste(loc.names[i], loc.results[i], sep="\t"))
  }
  cat("\nGlobal Estimate based on average heterozygosity\n")
  cat(summarise(x$global.het))
  if(!is.null(x$global.harm)){
    cat("\nGlobal Estimate based on harmonic mean of statistic\n")
    cat(summarise(x$global.harm))
  }
}

