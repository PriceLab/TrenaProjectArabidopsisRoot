library(RUnit)
source("./createTestMatrix.R")
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_createTestMatrix()

} # runTests
#----------------------------------------------------------------------------------------------------
test_createTestMatrix <- function()
{
   message(sprintf("--- test_createTestMatrix"))

} # test_createTestMatrix
#----------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()


