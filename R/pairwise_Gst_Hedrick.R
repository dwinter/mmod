#' Calculates pairwise values of Hedrick's G'st
#'
#' This function calculates Hedrick's G'st, a measure of genetic
#' differentiation, between all combinations of populaitons
#' in a genind object.
#'
#' @param x genind object (from package adegenet)
#' @param linearized logical, if TRUE will turned linearized G'st (1/()1-G'st))
#' @export
#' @examples
#' 
#' data(nancycats)
#' pairwise_Gst_Hedrick(nancycats[1:26,])
#' @references
#'  Hedrick, PW. (2005), A Standardized Genetic Differentiation Measure. Evolution 59: 1633-1638. 
#' @references
#'  Merimans, PG and Hedrick PW. (2010), Assessing population structure: FST and related measures. Molecular Ecology Resources 11: 5-18
#' @family pairwise
#' @family Hedrick

pairwise_Gst_Hedrick<- function(x, linearized=FALSE) {
  pops <- seppop(x)
  n.pops <- length(pops)
  #all combinations 
  allP <- combn(1:n.pops, 2)
  # calculate the statistic
  pair <- function(index.a,index.b){
    a <- pops[[index.a]]
    b <- pops[[index.b]]
    temp <- repool(a,b)
    return(Gst_Hedrick(temp)$global)
    }
  res <- sapply(1:dim(allP)[2], function(i) pair(allP[,i][1], allP[,i][2]))
  attributes(res) <- attributes(dist(1:n.pops))
  if(linearized){
   return(res/(1-res))
  }
  return(res)
}
