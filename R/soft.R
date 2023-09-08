soft <- function(x,d){
  return(sign(x)*pmax(0, abs(x)-d))
}
