##Sort the table by variance
tbl <- tbl.root
tbl$variance <- apply(tbl[, 3:10], 1, var)
tbl.sorted <- tbl[order(tbl$Alias, tbl$variance, decreasing= TRUE), ]

##Get rid of duplicates
tbl.sorted.nodupes <- tbl.sorted[-which(duplicated(tbl.sorted$Alias)), ]
rownames(tbl.sorted.nodupes) <- tbl.sorted.nodupes$Alias
tbl.sorted.nodupes <- tbl.sorted.nodupes[, -c(1:2)]

##Expand rows with multiple names
#----------------------------------------------------------------------------------------------------
library(RUnit)
#----------------------------------------------------------------------------------------------------
if(!exists("tbl.test")){
  print("fresh load of our data file")
  tbl.test <- tbl.sorted.nodupes[1:25, 1:10]
}
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
  test_splitMultipleNames()
  test_demo.doAllRows()
} # runTests
#----------------------------------------------------------------------------------------------------
splitMultipleNames <- function(tbl)
{
  tbl <- tbl.test[10,]
  stopifnot(nrow(tbl) == 1)    # we only handle data.frames with one row
  tokens <- strsplit(rownames(tbl), ";")[[1]]
  tbl.result <- tbl[rep(seq_len(nrow(tbl)), each=length(tokens)),]
  rownames(tbl.result) <- tokens
  
  return(tbl.result)
  
} # splitMultipleNames
#----------------------------------------------------------------------------------------------------
demo.doAllRows <- function(tbl)
{
  tbl <- tbl.test
  x <- lapply(1:nrow(tbl), function(i) splitMultipleNames(tbl[i,]))
  tbl.combined <- do.call(rbind, x)
  return(tbl.combined)
  
} # demo.doAllRows
#----------------------------------------------------------------------------------------------------
test_splitMultipleNames <- function()
{
  message(sprintf("--- test_splitMultipleNames"))   # announce what test you are running
  
  tbl.5 <- tbl.test[10,]
  tbl.fixed <- splitMultipleNames(tbl.5)
  checkEquals(nrow(tbl.fixed), 2)
  checkEquals(rownames(tbl.fixed), c("ATMG01250","AT2G07697"))
  
} #test_splitMultipleNames
#----------------------------------------------------------------------------------------------------
test_demo.doAllRows <- function()
{
  message(sprintf("---test_demo.doAllRows")) #announce what you're running
  tbl.expanded <- demo.doAllRows(tbl.test)
  
  checkTrue(nrow(tbl.expanded) > nrow(tbl.test))
  
  multiple.orf.names.orig <- length(grep(";", rownames(tbl.test)))
  checkEquals(multiple.orf.names.orig, 13)
  
  multiple.orf.names.expanded <- length(grep(";", rownames(tbl.expanded)))
  checkEquals(multiple.orf.names.expanded, 0)
  
} # test_splitMultiNamesToRows
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
  runTests()
