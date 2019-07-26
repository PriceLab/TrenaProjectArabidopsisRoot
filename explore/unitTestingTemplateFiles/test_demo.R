library(RUnit)
source("demo.R")
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_demo_1()
   test_demo_2()

} # runTests
#----------------------------------------------------------------------------------------------------
test_demo_1 <- function()
{
   message(sprintf("--- test_demo_1"))
   results <- demo(1, 2)
   checkEquals(results, 3)

} # test_demo_1
#----------------------------------------------------------------------------------------------------
test_demo_2 <- function()
{
   message(sprintf("--- test_demo_2"))
   results <- demo(-1, 1)
   checkEquals(results, 0)

} # test_demo_2
#----------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()


