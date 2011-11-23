#' Calculate differentiation stats for a jacknife sample of a Genind opject
#'
#' Calcutes
#'
#' @param x genind object (from package adegenet)
#' @param sample_frac fraction of pops to sample in each replication (default 0.5)
#' @param nreps number of jacknife replicates to run (default 1000)
#' @export
#' @examples
#'\dontrun{  
#' data(nancycats)
#' jacknife_populations(nancycats)
#' }

jacknife_populations <- function(x, sample_frac=0.5, nreps=1000){
 rep <- function(i,d){
  if( i %% 50 == 0){
     print(paste(i,"of", nreps,"reps completed"))
  }
  to.sample <- sample(d$pop.names, length(d$pop.names) * sample_frac)
  return(diff_stats(d[d$pop.names %in% to.sample,])$global)
 }
 return(sapply(1:nreps, rep, x))
}
  
