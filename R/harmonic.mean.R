harmonic.mean <- function(x){
  if(! all(x >= 0)){
    return(NA)   
    }
 return(1/mean(1/x))
}
