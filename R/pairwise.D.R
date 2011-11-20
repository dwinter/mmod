#' Calculates pairwise values of Jost's D
#'
#' This function calculates Jost's D, a measure of genetic
#' differentiation, between all combinations of populaitons
#' in genind object.
#'
#' @param x genind object (from package adegenet)
#' @export
#' @examples
#' 
#' data(nancycats)
#' pairwise.D(nancycats)

pairwise.D <- function(x) {
  pops <- seppop(x)
  n.pops <- length(pops)
  #all combinations 
  allP <- combn(1:n.pops, 2)
  # calculate tfh statistic
  pair <- function(index.a,index.b){
    a <- pops[[index.a]]
    b <- pops[[index.b]]
    temp <- repool(a,b)
    return(D.Jost(temp)$global.het)
    }
  res <- sapply(1:dim(allP)[2], function(i) pair(allP[,i][1], allP[,i][2]))
  attributes(res) <- attributes(dist(1:n.pops))
  return(as.matrix(res))
}
