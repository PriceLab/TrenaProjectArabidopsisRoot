library(RUnit)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_identifyMultiOrfNames()
   
} # runTests
#----------------------------------------------------------------------------------------------------
# create function identifyMultiOrfNames
identifyMultiOrfNames <- function(tbl)
{
  singleOrfNames <- unlist(strsplit(grep(";", tbl, value = TRUE), ";"))
  return(singleOrfNames)
}
#---------------------------------------------------------------------------------------------------
test_identifyMultiOrfNames <- function()
{
   message(sprintf("--- test_identifyMultiOrfNames"))
   
   result <- identifyMultiOrfNames("abcd;efgh")
   checkEquals(result, c("abcd", "efgh"))
   
   result.2 <- identifyMultiOrfNames("abcd;efgh;ijkl")
   checkEquals(result.2, c("abcd", "efgh", "ijkl"))
   
   result.3 <- identifyMultiOrfNames("abcd;efgh;ijkl;mnop")
   checkEquals(result.3, c("abcd", "efgh", "ijkl", "mnop"))
   
   result.4 <- identifyMultiOrfNames("abcd;efgh;ijkl;mnop;qrst")
   checkEquals(result.4, c("abcd", "efgh", "ijkl", "mnop", "qrst"))

} # test_identifyMultiOrfNames
#----------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
#---------------------------------------------------------------------------------------------------