Phi_st_meirmans <- function(x){


	add_max_distance <- function(x, index){
	sum(sapply(x[-1*index], function(y) y*x[index]))
	}	
	
	amova_st <- function(dm, dropped){
	
	if(is.null(dropped)){
		pops <- x@pop
	}
	else {
		pops <- droplevels(x@pop[-1 * dropped])
	}

	#values used multiple times
	sq_distance <- dm^2
	pop.freqs <- table(pops)
	n <- dim(dm)[1]
	df <- length(unique(pops)) - 1 
	within_dists <- sapply(unique(pops), function(p) sum(sq_distance[outer(pops==p, pops==p) == 1] ))
	ncoef <- (n - sum(pop.freqs^2)/n)/df

	
	#normal old AMOVA
	SSD_total <- sum(sq_distance/(2 * n))
	SSD_within <- sum(within_dists / (as.numeric(pop.freqs) *2))
	SSD_among <- SSD_total - SSD_within
  
	MSD_total <-SSD_total / (n-1)
	MSD_among <- SSD_among / df
	MSD_within <- SSD_within / ((n -1 ) - df)
	sigma2a <- (MSD_among - MSD_within)/ncoef
	phi <- sigma2a /(sigma2a + MSD_within)
	
	#Now AMOVA to the max
	max_between <- sapply(unique(pops), max_between_dist, pop.freqs)
	max_SSD_total <- sum(within_dists + max_between) / (2*n)
	max_MSD_among <- (max_SSD_total - SSD_within) / df
	sigma2_prime <- (max_MSD_among - MSD_within)/ncoef
	phi_max = sigma2a_prime /(sigma2a_prime + MSD_within)
	phi_prime <- phi/phi_max
	
	return(phi_prime)
	}

	global <- with(dist.codom(x), amova_st(distances, dropped))
	loci <- sapply(dist.codom(x, global=FALSE), function(d) with(d, amova_st(distances, dropped)))
	return(list(per.locus=loci, global=global))
}

