#' Calculates pairwise values of Nei's Gst
#'
#' This function calculates Nei's Gst, a measure of genetic
#' differentiation, between all combinations of populaitons
#' in a genind object.
#'
#' @return dist A distance matrix with between-population values of Gst
#' @param x genind object (from package adegenet)
#' @param linearized logical, if TRUE will turned linearized Gst (1/(1-Gst))
#' @export
#' @examples
#' 
#' data(nancycats)
#' pairwise_Gst_Nei(nancycats[1:26,])
#' @references
#'  Nei M. (1973) Analysis of gene diversity in subdivided populations. PNAS: 3321-3323. 
#' @references
#'  Nei M, Chesser RK. (1983). Estimation of fixation indices and gene diversities. Annals of Human Genetics. 47: 253-259.
#' @family pairwise
#' @family Nei


pairwise_Gst_Nei <- function(x, linearized=FALSE) {
  pops <- seppop(x)
  n.pops <- length(pops)
  #all combinations 
  allP <- utils::combn(1:n.pops, 2)
  # calculate tfh statistic
  pair <- function(index.a,index.b){
    a <- pops[[index.a]]
    b <- pops[[index.b]]
    temp <- repool(a,b)
    return(Gst_Nei(temp)$global)
    }
  res <- sapply(1:dim(allP)[2], function(i) pair(allP[,i][1], allP[,i][2]))
  attributes(res) <- list(class="dist", Diag=FALSE, Upper=FALSE, 
                          Labels=x@pop.names,Size=n.pops)
  if(linearized){
     return(res/(1-res))
  }
  return(res)
}
