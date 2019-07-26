library(RUnit)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_splitNames()
   
} # runTests
#----------------------------------------------------------------------------------------------------
# create function splitNames
splitNames <- function(string)
{
   if(!grepl(";", string))
      return(string)
  
  singleOrfNames <- unlist(strsplit(string, ";"))
  
  return(singleOrfNames)
}
#---------------------------------------------------------------------------------------------------
test_splitNames <- function()
{
   message(sprintf("--- test_splitNames"))
   checkEquals(splitNames("abc"), "abc")   
   checkEquals(splitNames(""), "")

   result <- splitNames("abcd;efgh")
   checkEquals(result, c("abcd", "efgh"))
   
   result.2 <- splitNames("abcd;efgh;ijkl")
   checkEquals(result.2, c("abcd", "efgh", "ijkl"))
   
   result.3 <- splitNames("abcd;efgh;ijkl;mnop")
   checkEquals(result.3, c("abcd", "efgh", "ijkl", "mnop"))
   
   result.4 <- splitNames("abcd;efgh;ijkl;mnop;qrst")
   checkEquals(result.4, c("abcd", "efgh", "ijkl", "mnop", "qrst"))

} # test_splitNames
#----------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
#---------------------------------------------------------------------------------------------------

