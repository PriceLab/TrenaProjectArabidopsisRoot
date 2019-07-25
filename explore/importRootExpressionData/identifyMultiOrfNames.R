library(RUnit)
#----------------------------------------------------------------------------------------------------
# create function identifyMultiOrfNames
identifyMultiOrfNames <- function(tbl)
{
  tbl <- read.table("roots.tsv", sep = "\t", header= TRUE, as.is=TRUE, nrow=400)
  multiOrfNames <- grep(";", tbl$Alias, value =TRUE)
  singleOrfNames <- unlist(strsplit(multiOrfNames, ";"))
}  
