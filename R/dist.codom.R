#' Calculate distance between individual for co-dominant locus
#'
#' This function calculates the distance between individuals in a genind
#' object based on their genotypes. 
#' Specifically, the simple metric of Kosman and Leonard (2005) in which 
#' distance is calculated as a propotion of shared alleles at each locus.
#' 
#'
#' @param x genind object (from package adegenet)
#' @param matrix boolean: if TRUE return matrix (dist object if FALSE)
#' @param global boolean: if TRUE, return a single global estimate based on all
#' loci. If FALSE return a list of matrices for each locus.
#' if FALSE 
#' @param na.rm boolean: if TRUE remove individuals with NAs
#'
#' @return either a list of distance matrices, one for each locus or a single 
#' matrix containing the mean distance between individuals across all loci
#' @return Dropped for each distance matrix and object of class "na.action" 
#' containing indices to those indivudals in the genind object which where 
#' omitted due to having NAs
#' @export
#' 
#' @references
#'  Kosman E., Leonard, K.J. Similarity coefficients for molecular markers in studies of genetic relationships between individuals for haploid
#'  diploid, and polyploid species. Molecular Ecology. 14: 415-424
#'  
#' @examples
#' data(nancycats)
#' dist.codom(nancycats[,1])


dist.codom <- function(x, matrix=TRUE, global=TRUE, na.rm=TRUE){
  
  if (! ploidy(x) == 2){
   stop("Sorry, the function dist.codom only works diploid datasets at the moment. Fixing this is on the TODO list, contact the author if this causes a problem...")
  }

  pair <- function(l, x, y){
  #find distance between individuls x and y for locus l
    return( sum( abs(l@tab[x,]- l@tab[y,]) ) / 2 )
  }
  
  per.loc <- function(l){
	if(na.rm){ 
	  l@tab <- na.omit(l@tab)
	  dropped <- attr(l@tab, "na.action") 
	}
	n.inds <- dim(l@tab)[1]
	allP <- combn(1:n.inds, 2)
	res <- sapply(1:dim(allP)[2], function(i) pair(l, allP[,i][1], allP[,i][2]))
	#trick to turn vector into distance matrix
	attributes(res) <- attributes(dist(1:n.inds))	
	if(matrix) {
		res <- as.matrix(res)
		}
	if(na.rm){
		return(list(distances=res, dropped=dropped))
		}
	else return(res)		
	}

  
	if(global){
		has_nas <- FALSE
		if(any(is.na(x@tab))){
		# if we want a global estimate make new genind with out NAs
		# (just chopping out NAs from x@tab won't work)
			has_nas  <- TRUE
			warning("removing all individuals with NAs for global estimate")
			n_tab <- na.omit(x@tab)
			dropped <- attr(n_tab, "na.action")
			n_pop <- droplevels(x@pop[ -1 * dropped])
			x <- genind(tab=n_tab, pop=n_pop)
		}	
		dists <- lapply(seploc(x), per.loc)
		xbar <- Reduce("+", lapply(dists, function(x) x$dist))  / length(dists)
		if(has_nas){
			return(list(distances=xbar, dropped=dropped))
		}
		else return(list(distances=xbar,dropped=NULL))
	}
	
	else{
		dists <- lapply(seploc(x), per.loc)	
		return(dists)
	}
}
