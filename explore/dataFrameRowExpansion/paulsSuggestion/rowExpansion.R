library(RUnit)
#----------------------------------------------------------------------------------------------------
if(!exists("tbl.test")){
    tbl.test <- tbl.sorted.nodupe.10by10
    }
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
  test_splitMultipleNamesIntoTheirOwnRows()

} # runTests
#----------------------------------------------------------------------------------------------------
tbl <- tbl.sorted.nodupe.10by10[10,]
df <- tbl 
splitMultipleNamesIntoTheirOwnRows <- function(tbl)
{
  stopifnot(nrow(tbl) == 1)    # we only handle data.frames with one row
  
  if(!grepl(";", rownames(tbl))){
     return(tbl)  # just give it back to the caller
     }

  tokens <- strsplit(rownames(tbl), ";")[[1]]
  df[rep(seq_len(nrow(df)), each=length(tokens)),]
  
  
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

} # splitMultipleNamesIntoTheirOwnRows
#----------------------------------------------------------------------------------------------------
test_splitMultipleNamesIntoTheirOwnRows <- function()
{
  message(sprintf("--- test_splitMultipleNamesIntoTheirOwnRows"))   # announce what test you are running

     #------------------------------------------------------------------------
     #  do the simplest test first: a row with a simple name (no semicolons)
     #------------------------------------------------------------------------

  tbl.1 <- tbl.test[1,]
  tbl.fixed <- splitMultipleNamesIntoTheirOwnRows(tbl.1)
  checkEquals(nrow(tbl.fixed), 1)
  checkEquals(rownames(tbl.fixed), rownames(tbl.1))

  tbl.2 <- tbl.test[10,]
  tbl.fixed <- splitMultipleNamesIntoTheirOwnRows(tbl.2)
  checkEquals(nrow(tbl.fixed), 2)
  checkEquals(rownames(tbl.fixed),c("ATMG01250", "AT2G07697"))
  
  tbl.3 <- tbl.1
  rownames(tbl.3) <- "aaa;bbb;ccc"
  tbl.fixed <- splitMultipleNamesIntoTheirOwnRows(tbl.3)
  checkEquals(nrow(tbl.fixed), 3)
  checkEquals(rownames(tbl.fixed),c("aaa", "bbb", "ccc"))
  
  tbl.4 <- tbl.1
  rownames(tbl.4) <- "abc;def;ghi;jkl"
  tbl.fixed <- splitMultipleNamesIntoTheirOwnRows(tbl.4)
  checkEquals(nrow(tbl.fixed), 4)
  checkEquals(rownames(tbl.fixed),c("abc", "def", "ghi", "jkl"))
  

  
} # test_splitMultiple
#----------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
