#' Produce bootstrap samples from each subpopulation of a genind object
#'
#' This function produces bootstrap samples from a genind object, with each
#' subpopulation resampled according to its size. Because there are many 
#' statistics that you may wish to calculte from these samples, this function
#' returns a list of genind objects representing bootsrap samples that can then
#' be futher processed (see examples).
#' 
#'
#' @param x genind object (from package adegenet)
#' @param nreps numeric number of bootstrap replicates to perform (default 1000)
#' @export
#' @examples
#'\dontrun{  
#' data(nancycats)
#' bs <- chao_bootstrap(nancycats)
#' bs_D <- sapply(bs, D_Jost)
#' hist(unlist(bs_D[2,]))
#' quantile( unlist(bs_D[2,]), c(0.025, 0.5, 0.975) )
#'}
#' @references
#' 
#'
  
chao_bootstrap <- function(x, nreps = 1000){
  one_rep <- function(){
    temp <- x
    bs <- by(x@tab, pop(x),function(x) x[sample(1:dim(x)[1],replace=T),] )
    temp@tab <- as.matrix(do.call("rbind", bs))
    return(temp)
  }
  res <- replicate(nreps, one_rep() )
  return(res)
}


summarise_bootsrap(bs, statistic, round_to=4, plot=FALSE){
  nreps <- dim(bs)[2]
  res <- sapply(bs, statistic)
  per.locus <- matrix(unlist(res[1,]), nrow=nreps)
  global.het <- unlist(res[2,])
  #Only D_Jost has another estimate of global value
  if(identical(stat, D_Jost)){
    global.harm <- unlist(res[3,])
  }
  
  summarise <- function(x){
    ci <- round(quantile(x, c(0.025, 0.975)), round_to)
    mean <- round(mean(x), round_to)
    line <- paste(mean, " (", ci[1], "-", ci[2], ")", sep="")  
    return(line)
  }
  apply(per.locus, 2, summarise)
  sapply(global.het, summarise)

}
