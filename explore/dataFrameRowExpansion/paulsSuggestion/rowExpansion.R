library(RUnit)
#----------------------------------------------------------------------------------------------------
if(!exists("tbl.test")){
  print("fresh load of our data file")
    tbl.test <- get(load(file= "10by10rootTable.RData"))
    }
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
  test_splitMultipleNamesIntoTheirOwnRows()

} # runTests
#----------------------------------------------------------------------------------------------------
obsolete.splitMultipleNamesIntoTheirOwnRows <- function(tbl)
{
  stopifnot(nrow(tbl) == 1)    # we only handle data.frames with one row
  
  if(!grepl(";", rownames(tbl))){
     return(tbl)  # just give it back to the caller
     }

  tokens <- strsplit(rownames(tbl), ";")[[1]]
  tbl[rep(seq_len(nrow(tbl)), each=length(tokens)),]
  
  
  if(length(tokens) == 2)
     tbl.result <- rbind(tbl, tbl)

  if(length(tokens) == 3)
     tbl.result <- rbind(tbl, tbl, tbl)

  if(length(tokens) == 4)
     tbl.result <- rbind(tbl, tbl, tbl, tbl)

  rownames(tbl.result) <- tokens

  return(tbl.result)
  # your code to split the rowname, and duplicate the row values, building
  # up new data.frame, of two or more rows,  will happen here

} # obsolete.splitMultipleNamesIntoTheirOwnRows
#----------------------------------------------------------------------------------------------------
splitMultipleNames <- function(tbl)
{
  stopifnot(nrow(tbl) == 1)    # we only handle data.frames with one row
  tokens <- strsplit(rownames(tbl), ";")[[1]]
  tbl.result <- tbl[rep(seq_len(nrow(tbl)), each=length(tokens)),]
  rownames(tbl.result) <- tokens
  
  return(tbl.result)

} # splitMultipleNames
#----------------------------------------------------------------------------------------------------
test_splitMultipleNames <- function()
{
  message(sprintf("--- test_splitMultipleNamesIntoTheirOwnRows"))   # announce what test you are running

     #------------------------------------------------------------------------
     #  do the simplest test first: a row with a simple name (no semicolons)
     #------------------------------------------------------------------------

  tbl.1 <- tbl.test[1,]
  tbl.fixed <- splitMultipleNames(tbl.1)
  checkEquals(nrow(tbl.fixed), 1)
  checkEquals(rownames(tbl.fixed), rownames(tbl.1))

  tbl.2 <- tbl.test[10,]
  tbl.fixed <- splitMultipleNames(tbl.2)
  checkEquals(nrow(tbl.fixed), 2)
  checkEquals(rownames(tbl.fixed),c("ATMG01250", "AT2G07697"))
  
  tbl.3 <- tbl.1
  rownames(tbl.3) <- "aaa;bbb;ccc"
  tbl.fixed <- splitMultipleNames(tbl.3)
  checkEquals(nrow(tbl.fixed), 3)
  checkEquals(rownames(tbl.fixed),c("aaa", "bbb", "ccc"))
  
  tbl.4 <- tbl.1
  rownames(tbl.4) <- "abc;def;ghi;jkl"
  tbl.fixed <- splitMultipleNames(tbl.4)
  checkEquals(nrow(tbl.fixed), 4)
  checkEquals(rownames(tbl.fixed),c("abc", "def", "ghi", "jkl"))
  
} # test_splitMultiple
#----------------------------------------------------------------------------------------------------
test_MultiRowTableSplitMultipleNames <- function() 
{  
  message(sprintf("--- test_splitMultiple"))   # announce what test you are running
  
  tbl.5 <- tbl.test[8,]
  tbl.fixed <- splitMultipleNames(tbl.5)
  checkEquals(nrow(tbl.fixed), 5)
  checkEquals(rownames(tbl.fixed), rownames(tbl.5))
} 
#----------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
