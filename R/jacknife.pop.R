jacknife.populations <- function(x, sample_frac=0.5, nreps=1000){
 rep <- function(i,d){
  if( i %% 50 == 0){
     print(paste(i,"of", nreps,"reps completed"))
  }
  to.sample <- sample(d$pop.names, length(d$pop.names) * sample_frac)
  return(diff.stats(d[d$pop.names %in% to.sample,])$global)
 }
 return(sapply(1:nreps, rep, x))
}
  
