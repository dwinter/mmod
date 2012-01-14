#' Harmonic mean
#'
#' Calculate the harmonic mean of a numeric vector
#' (will return NA if there are any negative numbers in the vector)
#' @param x numeric vector
#' @export
#' @examples
#' 
#' data(nancycats)
#' pop.sizes <- table(pop(nancycats))
#' harmonic_mean(pop.sizes)


harmonic_mean <- function(x){
  if(! all(x >= 0)){
    return(NA)   
    }
 return(1/mean(1/x))
}
