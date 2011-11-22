#' Harmonic mean
#'
#' Calculate the harmonic mean of a numeric vector
#'
#' @param x numeric vector
#' @export
#' @examples
#' 
#' data(nancycats)
#' pop.sizes <- table(pop(nancycats))
#' harmonic.mean(pop.sizes)


harmonic.mean <- function(x){
  if(! all(x >= 0)){
    return(NA)   
    }
 return(1/mean(1/x))
}
