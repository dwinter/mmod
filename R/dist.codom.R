#' Calculate distance between individual for co-dominant locus
#'
#' This function calculates the distance between individuals in a genind
#' object based on their genotypes. 
#' Specifically, the simple metric of Kosman and Leonard (2005) in which 
#' distance is calculated as a propotion of shared alleles at each locus.
#' 
#'
#' @param x genind object (from package adegenet)
#'
#' @return either a list of distance matrices, one for each locus or a single 
#' matrix containing the mean distance between individuals across all loci
#' @export
#' 
#' @references
#'  Kosman E., Leonard, K.J. Similarity coefficients for molecular markers in studies of genetic relationships between individuals for haploid
#'  diploid, and polyploid species. Molecular Ecology. 14: 415-424
#'  
#' @examples
#' data(nancycats)
#' dist_codom(nancycats[,1])

dist.codom <- function(g, matrix=TRUE, seplocs=TRUE){
  per.locus <- function(l){
      n.inds <- dim(l@tab)[1]
      allP <- combn(1:n.inds, 2)
      pair <- function(x, y){
        return( sum( abs(l@tab[x,]- l@tab[y,]) ) / 2 )
      }
      res <- sapply(1:dim(allP)[2], function(i) pair(allP[,i][1], allP[,i][2]))
      #trick to turn vector into distance matrix
      attributes(res) <- attributes(dist(1:n.inds))
      if(!matrix){
	    return(res)
      }
      else return(as.matrix(res))
  }
  loci <- lapply(seploc(g), per.locus)  
  if (!seplocs){
    if(any(is.na(g@tab))){
        warning("Ommitting individuals with NAs")
        no_nas <- apply(g@tab, 1, function(x) !any(is.na(x)) )
        loci <- lapply(loci, function(x) x[no_nas,no_nas])
    }
    return( Reduce("+", loci) / length(loci) )
  }  
  else return(loci)
}
