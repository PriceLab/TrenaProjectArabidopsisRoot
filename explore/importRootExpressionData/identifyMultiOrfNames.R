library(RUnit)
#----------------------------------------------------------------------------------------------------
# create function identifyMultiOrfNames
identifyMultiOrfNames <- function(tbl)
{
  tbl <- read.table("roots.tsv", sep = "\t", header= TRUE, as.is=TRUE, nrow=400)
  multiOrfNames <- grep(";", tbl$Alias, value =TRUE)
  singleOrfNames <- unlist(strsplit(multiOrfNames, ";"))
}
#---------------------------------------------------------------------------------------------------
test_identifyMultiOrfNames <- function()
{
   message(sprintf("--- test_identifyMultiOrfNames"))
   result <- idenitfyMultiOrfNames("abcd;efgh")
   checkEquals(result, c("abcd", "efgh"))


} # test_identifyMultiOrfNames
#---------------------------------------------------------------------------------------------------
