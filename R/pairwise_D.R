#' Calculates pairwise values of Jost's D
#'
#' This function calculates Jost's D, a measure of genetic
#' differentiation, between all combinations of populaitons
#' in a genind object.
#'
#' @param x genind object (from package adegenet)
#' @export
#' @examples
#' 
#' data(nancycats)
#' pairwise_D(nancycats[1:26,])



fast_repool <- function(...){
    x <- list(...)
    listTab <- lapply(x, genind2df, usepop = FALSE)
    getPop <- function(obj) {
        pop <- obj$pop
        levels(pop) <- obj$pop.names
        return(pop)
    }
    listPop <- lapply(x, getPop)
    pop <- unlist(listPop, use.name = FALSE)
    pop <- factor(pop)
    markNames <- colnames(listTab[[1]])
    listTab <- lapply(listTab, function(tab) tab[, markNames, 
        drop = FALSE])
    tab <- listTab[[1]]
    for (i in 2:length(x)) {
        tab <- rbind(tab, listTab[[i]])
    }
    res <- df2genind(tab, pop = pop, ploidy = x[[1]]@ploidy, 
        type = x[[1]]@type)
    res$call <- match.call()
    return(res)
}


pairwise_D <- function(x) {
  pops <- seppop(x)
  n.pops <- length(pops)
  #all combinations 
  allP <- combn(1:n.pops, 2)
  # calculate tfh statistic
  pair <- function(index.a,index.b){
    a <- pops[[index.a]]
    b <- pops[[index.b]]
    temp <- fast_repool(a,b)
    return(D_Jost(temp)$global.het)
    }
  res <- sapply(1:dim(allP)[2], function(i) pair(allP[,i][1], allP[,i][2]))
  attributes(res) <- attributes(dist(1:n.pops))
  return(as.matrix(res))
}





