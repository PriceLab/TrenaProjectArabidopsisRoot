library(RUnit)
#----------------------------------------------------------------------------------------------------
addFour.old <- function(d)
{
  stopifnot(is.numeric(d))

  sum <- d + 4

  return(sum)

} # addFour.old
#----------------------------------------------------------------------------------------------------
addFour <- function(a)
{
   return(add(a, 4))

} # addFour
#----------------------------------------------------------------------------------------------------
add <- function(a, b)
{
  stopifnot(is.numeric(a))
  stopifnot(is.numeric(b))

  return(a + b)

} # add
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
  test_addFour()

} # runTests
#----------------------------------------------------------------------------------------------------
test_addFour <- function()
{
  message(sprintf("--- test_addFour"))

  checkException(addFour("doofus"), silent=TRUE)

  checkEquals(addFour(28), 32)
  #checkEquals(addFour.v2(28), 32)

  checkEquals(addFour(-1000), -996)
  checkException(addFour("8"), silent=TRUE)

} # test_createTestMatrix
#----------------------------------------------------------------------------------------------------
if(!interactive())
  runTests()

