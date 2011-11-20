pairwise.Hedrick.Nei <- function(x) {
  pops <- seppop(x)
  n.pops <- length(pops)
  #all combinations 
  allP <- combn(1:n.pops, 2)
  # calculate tfh statistic
  pair <- function(index.a,index.b){
    a <- pops[[index.a]]
    b <- pops[[index.b]]
    temp <- repool(a,b)
    return(Gst.Hedrick(temp)$global)
    }
  res <- sapply(1:dim(allP)[2], function(i) pair(allP[,i][1], allP[,i][2]))
  attributes(res) <- attributes(dist(1:n.pops))
  return(as.matrix(res))
}
