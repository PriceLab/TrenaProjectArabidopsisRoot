library(RUnit)

addFour <- function(d)
{
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
test_addFour <- function(d)
{
  message(sprintf("--- test_addFour"))
  
  checkTrue(addFour(28) > 28)
  checkEquals(addFour(28), 32)
  checkTrue(class(addFour(28))== "numeric")
  
} # test_createTestMatrix
#----------------------------------------------------------------------------------------------------
if(!interactive())
  runTests()

