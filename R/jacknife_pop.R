#' Calculate differentiation stats for a jacknife sample of a Genind opject
#'
#' Makes a series of jacknife samples across populations from a Genind object 
#' and calculates differentiation stats for each sample.
#'
#' @param x genind object (from package adegenet)
#' @param sample_frac fraction of pops to sample in each replication (default 0.5)
#' @param nreps number of jacknife replicates to run (default 1000)
#' @export
#' @examples
#'\dontrun{  
#' data(nancycats)
#' obs <- diff_stats(nancycats)
#' jn <- jacknife_populations(nancycats)
#' D_sampled <-jn[5,]
#' hist(D_sampled)
#' abline(h=obs$global)
#' }

jacknife_populations <- function(x, sample_frac=0.5, nreps=1000){
 rep <- function(i,d){
  if( i %% 50 == 0){
     cat(paste(i,"of", nreps,"reps completed"))
  }
  to.sample <- sample(d$pop.names, length(d$pop.names) * sample_frac)
  return(diff_stats(d[d$pop.names %in% to.sample,])$global)
 }
 return(sapply(1:nreps, rep, x))
}
  
