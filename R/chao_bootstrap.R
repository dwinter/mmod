#' Produce bootstrap samples from each subpopulation of a genind object
#'
#' This function produces bootstrap samples from a genind object, with each
#' subpopulation resampled according to its size. Because there are many 
#' statistics that you may wish to calculte from these samples, this function
#' returns a list that can be futher processed (see examples).
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


