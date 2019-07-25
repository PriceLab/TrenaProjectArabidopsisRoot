addFour <- function(d)
{
  d <- 3
  stopifnot(class(d)== "numeric")
  sum <- d + 4
  return(sum)
  
} #addFour()