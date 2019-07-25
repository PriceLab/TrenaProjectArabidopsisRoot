addFour <- function(d)
{
  d <- 3
  stopifnot(class(d)== "numeric")
  sum <- d + 4
  return(sum)
  
} #addFour()

#----------------------------------------------------------------------------------------------------
runTests <- function()
{
  test_addFour()
  
} # runTests
#----------------------------------------------------------------------------------------------------
test_addFour <- function()
{
  message(sprintf("--- test_addFour"))
  
  sum <- addFour(d)
  checkTrue(sum > d)
  checkEquals(sum, d + 4)
  checkTrue(class(sum)== "numeric")
  
} # test_createTestMatrix
#----------------------------------------------------------------------------------------------------
if(!interactive())
  runTests()
