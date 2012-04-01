#Internal function, used to get maximum possible distance for dataset
max_between_dist <- function(pop_name, pop_freqs){
    inter_dists <- sapply(pop_freqs[names(pop_freqs) != pop_name], 
                      function(x) x * pop_freqs[pop_name])
    return(sum(inter_dists))
}

Phi_st_Meirmans <- function(x){
    
    amova_st <- function(dm){    
    
      if(any(is.na(locus@tab))){
        warning(paste("Ommitting NAs from locus", locus@loc.names))
        locus@tab <- na.omit(locus@tab)
        locus@pop <- droplevels(locus@pop[-1 * attr(locus@tab, "na.action")])
      }
      
      dm <- dist.codom(locus)^2
      n <- dim(dm)[1]
      pop.freqs <- table(pops)
      
      #values we will need multiple times  
      df <- length(unique(pops)) - 1
      # within populaton distances for each population
      # more readble to do this with melt/reshape of l@tab but 
      # avoid dependency
      within_dists <- sapply(unique(pops), function(p) 
                      sum(dm[outer(pops==p, pops==p) == 1] ))
      ncoef <- (n - sum(pop.freqs^2)/n)/df

      #Normal old AMOVA
      SSD_total <- sum(dm/(2 * n))
      SSD_within <- sum(within_dists / (as.numeric(pop.freqs) * 2))
      SSD_among <- SSD_total - SSD_within
      
      MSD_total <-SSD_total / (n-1)
      MSD_among <- SSD_among / df
      MSD_within <- SSD_within / ((n -1 ) - df)
      sigma2 <- (MSD_among - MSD_within)/ncoef
      phi = sigma2/(sigma2  + MSD_within)
      
      #Now compared to maximum value of SSD
      max_between <- sapply(unique(pops), max_between_dist, pop.freqs)
      max_SSD_total <- sum(within_dists + max_between) / (2*n)
      SSD_among_prime <- max_SSD_total - SSD_within
      MSD_among_prime <- SSD_among_prime / df
      sigma2_prime <- (MSD_among_prime - MSD_within)/ncoef
      phi_max = sigma2_prime/(sigma2_prime  + MSD_within)
      phi_prime = phi/phi_max
      
      return(c(SSD=SSD_among, SSD_prime = SSD_among_prime, df=df, 
               phi_st=phi, phi_prime_st=phi_prime))
     }
     phi_loci <- sapply(seploc(x), amova_st))
     phi_global <- mean(phi_loci[,5])
     return(loci=phi_loci, global=phi_global)
     
 }
 
 
